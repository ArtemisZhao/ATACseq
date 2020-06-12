library(BiocManager)
#BiocManager::install("ChIPQC")
library(ChIPQC)
#BiocManager::install("ChIPseeker")
library(ChIPseeker)
#BiocManager::install("limma")
library(limma)
library(DT)
library(tidyr)
library(Rsubread)
library(dplyr)
library(rtracklayer)
#BiocManager::install("DESeq2")
library(DESeq2)
library(soGGi)
library(Rsamtools)
library(GenomicAlignments)
library(latticeExtra)
library(argparser,quietly = TRUE)
library(diffloop)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mygene)
library(tracktables)

p<-arg_parser("Differential expression analysis - DEseq")
p<-add_argument(p,"peak_type", help="Table contains peak type names for each file")
p<-add_argument(p, "group_name",help="Table contains covariate factor")
p<-add_argument(p, "blkList", help="black list bed file")
p<-add_argument(p, "--peak_dir", short="-pdir", default=".", help='Peak file directory')
p<-add_argument(p,"--bam_dir",short="-bdir", default=".", help='Bam file directory')
p<-add_argument(p, "--output_dir",short="-o",default=".",help="Output directory")
p<-add_argument(p, "--add_cov", default=NULL, help="Additional covariate, e.g. signal to noise ratio")
argv<-parse_args(p)


cat("Identifying a set of non-redudant peaks ...\n")

peaks <- dir(argv$peak_dir, pattern = "*.narrowPeak", full.names = TRUE)

myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

peakname<- t(read.table(argv$peak_type))
print(length(peakname))
print(length(myPeaks))

if (length(peakname)!=length(myPeaks)){
  print("check your name file")}

if(length(peakname)==length(myPeaks)) {
  names(myPeaks) <- peakname
}

Group <- factor(t(read.table(argv$group_name))) 

if (!argv$add_cov=="NULL"){
AddGroup <- factor(t(read.table(argv$add_cov)))
}

#consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")
myGRangesList<-GRangesList(myPeaks)   
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))

consensusToCount<-reducedConsensus


cat("Blacklisted region removal ... \n")

blklist <- import.bed(argv$blkList)

###################################### this step add chr to seqnames
consensusToCount<- addchr(consensusToCount)
consensusToCount <- consensusToCount[!consensusToCount %over% blklist ]

save(consensusToCount,file=paste0(argv$output_dir,"/consensusToCount.Rdata",sep=""))

cat("Keeping peaks present in more than two replicates ...\n ")

occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>%
  rowSums

#table(occurrences) %>% rev %>% cumsum
consensusToCount <- consensusToCount[occurrences >= 2, ] 

save(consensusToCount,file=paste0(argv$output_dir,"/consensusToCount2.Rdata",sep=""))

cat("Loading in bam files ...\n")
bamsToCount <- dir(argv$bam_dir, full.names = TRUE, pattern = "*.\\.bam$")

# indexBam(bamsToCount)

regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount),
              start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount),
              Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, 
                           isPairedEnd = TRUE, countMultiMappingReads = FALSE, maxFragLength = 100)

myCounts <- fcResults$counts

colnames(myCounts) <- peakname

save(myCounts, file = paste0(argv$output_dir,"/countsFromATAC.RData"))

cat("DEseq differential analysis ...\n")
if (!argv$add_cov=="NULL"){
  metaData <- data.frame(cbind(Group,AddGroup), row.names = colnames(myCounts))
  atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group+AddGroup+Group:AddGroup, 
                                    rowRanges = consensusToCount) 
  atacDDS <- DESeq(atacDDS)
  atac_Rlog <- rlog(atacDDS)
  
  save(atacDDS,file=paste0(argv$output_dir,"/atacDDS.Rdata"))
  
  } else{
  metaData <- data.frame(cbind(Group), row.names = colnames(myCounts))
 atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount) 
 atacDDS <- DESeq(atacDDS)
 atac_Rlog <- rlog(atacDDS)

 save(atacDDS,file=paste0(argv$output_dir,"/atacDDS.Rdata"))
}

resultname1<-resultsNames(atacDDS)

resultname<- gsub("vs_","",resultname)
resultname <- strsplit(resultname,c("_"))[-1]


a <- lapply(resultname[1:2],
            function(x){
              res <- results(atacDDS, contrast=x, format = "GRanges") 
              res <- res[order(res$pvalue)]
              res <- res[(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]
              gr <- annotatePeak(res,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
              gr<-as.data.frame(as.GRanges(gr))
              GeneName <- queryMany(gr$geneId, scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
              gr$gene.symbol <- GeneName$response$symbol
              gr
            })

for (i in 1:length(a)){
  write.table(a[[i]], file=paste0(argv$output_dir,"/",resultname1[i+1],".txt"), quote=F, sep="\t", row.names=T, col.names=T)
}
