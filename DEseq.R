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

p<-arg_parser("Differential expression analysis")
p<-add_argument(p,"peak_type", help="Table contains peak type names")
p<-add_argument(p, "group_name",help="Table contains group names")
p<-add_argument(p, "blkList", help="black list bed file")
p<-add_argument(p, "--peak_dir", short="-pdir", default=".", help='Peak file directory')
p<-add_argument(p,"--bam_dir",short="-bdir", default=".", help='Bam file directory')
p<-add_argument(p, "--output_dir",short="-o",default=".",help="Output directory")
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

metaData <- data.frame(Group, row.names = colnames(myCounts))


cat("DEseq differential analysis ...\n")
atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount) 
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)

save(atacDDS,file=paste0(argv$output_dir,"/atacDDS.Rdata"))

#pdf(file=paste0(argv$output_dir,"/PCA.pdf",width=6,height=6,bg="white",point=16))
#plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))
#dev.off()

###group need to specify
#AffMinusUnaff <- results(atacDDS, c("Group", "0h","24h"), format = "GRanges") 
#AffMinusUnaff <- AffMinusUnaff[order(AffMinusUnaff$pvalue)]

#print(AffMinusUnaff)
#save(AffMinusUnaff,file="pvalues.Rdata")

