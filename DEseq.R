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
p<-add_argument(p,"--cont",help="Whether the additional covariate is continuous",flag=TRUE)
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

Group <- read.table(argv$group_name,header=T)

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

#save(consensusToCount,file=paste0(argv$output_dir,"/consensusToCount.Rdata",sep=""))

cat("Keeping peaks present in more than two replicates ...\n ")

occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>%
  rowSums

#table(occurrences) %>% rev %>% cumsum
consensusToCount <- consensusToCount[occurrences >= 2, ] 

save(consensusToCount,file=paste0(argv$output_dir,"/consensusToCount.Rdata",sep=""))

cat("Loading in bam files ...\n")
bamsToCount <- dir(argv$bam_dir, full.names = TRUE, pattern = "*.\\.bam$")

# indexBam(bamsToCount)

cat("Calculating counts landing in peaks")
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount),
              start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount),
              Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, 
                           isPairedEnd = TRUE, countMultiMappingReads = FALSE, maxFragLength = 100)

myCounts <- fcResults$counts

colnames(myCounts) <- peakname

save(myCounts, file = paste0(argv$output_dir,"/countMatrix.RData"))


##############################################################################
 cat("DEseq differential analysis ...\n")

 ###contains a hard coding part
 atacDDS <- DESeqDataSetFromMatrix(myCounts, Group, ~Psphen+Celltype+
                                      Treat+skinhoming+Agebin+Celltype:Treat+
                                      Celltype:skinhoming+skinhoming:Treat+Gender, 
                                    rowRanges = consensusToCount) 
  atacDDS <- DESeq(atacDDS)
  atac_Rlog <- vst(atacDDS,blind=TRUE)
  
  save(atacDDS,file=paste0(argv$output_dir,"/atacDDS.Rdata"))
  

 #pdf(file=paste0(argv$output_dir,"atacPCA.pdf"),width=6,height=6.6,bg="white",pointsize=14)
 #plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))
 #dev.off()

 resultname1<-resultsNames(atacDDS)

 resultname <- resultname1[-1]


 a <- lapply(resultname,
            function(x){
              res <- results(atacDDS, name=x, format = "GRanges") 
              res <- res[order(res$pvalue)]
              res <- res[(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]
              gr <- annotatePeak(res,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
              gr<-as.data.frame(as.GRanges(gr))
              GeneName <- queryMany(gr$geneId, scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
              gr$gene.symbol <- GeneName$response$symbol
              gr
            })

for (i in 1:length(a)){
  write.table(a[[i]], file=paste0(argv$output_dir,"/",resultname1[i],".txt"), quote=F, sep="\t", row.names=T, col.names=T)
}

 
 
