library(ChIPQC)
library(ChIPseeker)
#BiocManager::install("limma")
library(limma)
library(Rsubread)
library(dplyr)
BiocManager::install("DESeq2")
library(DESeq2)


p<-arg_parser("Differential expression analysis")
p<-add_argument(p,"peak_type", help="Table contains peak type names")
p<-add_argument(p, "group_name",help="Table contains group names")
p<-add_argument(p, "blkList", help="black list bed file")
p<-add_argument(p, "--peak_dir", short="-pdir", default=".", help='Peak file directory')
p<-add_argument(p,"--bam_dir",short="-bdir", default=".", help='Bam file directory')
p<-add_argument(p, "--output_dir",short="-o",default=".",help="Output directory")
argv<-parse_args(p)


cat("Identifying a set of non-redudant peaks ...")

peaks <- dir(argv$dir, pattern = "*.narrowPeak", full.names = TRUE)

myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

peakname<- t(read.table(argv$peak_type))

###maywrong
if (length(peakname)!=length(myPeaks)){
  print("check your name file")}

if(length(peakname)==length(myPeaks)) {
  names(myPeaks) <- peaknames
}
Group <- factor(t(read.table(argv$group_name))) 

consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")

blklist <- import.bed(argv$blkList)

consensusToCount <- consensusToCount[!consensusToCount %over% blklist ]

save(consensusToCount,file=paste0(argv$output_dir,"/consensusToCount.Rdata",sep=""))

cat("Keeping peaks present in more than two replicates ")
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>%
  rowSums
table(occurrences) %>% rev %>% cumsum

consensusToCount <- consensusToCount[occurrences >= 2, ] 

bamsToCount <- dir(argv$bam_dir, full.names = TRUE, pattern = "*.\\.bam $")

# indexBam(bamsToCount)

regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount),
              start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount),
              Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, 
                           isPairedEnd = TRUE, countMultiMappingReads = FALSE, maxFragLength = 100)

myCounts <- fcResults$counts

colnames(myCounts) <- peaknames

save(myCounts, file = paste0(argv$output_dir,"/countsFromATAC.RData"))

metaData <- data.frame(Group, row.names = colnames(myCounts))

atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount) 
atacDDS <- DESeq(atacDDS)
atac_Rlog <- rlog(atacDDS)

plotPCA(atac_Rlog, intgroup = "Group", ntop = nrow(atac_Rlog))

###need to change
AffMinusUnaff <- results(atacDDS, c("Group", "Affected", "Unaffected"), format = "GRanges") 
AffMinusUnaff <- AffMinusUnaff[order(AffMinusUnaff$pvalue)]
print(AffMinusUnaff)


