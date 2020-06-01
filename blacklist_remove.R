library(BiocManager,quietly = TRUE)
#BiocManager::install("ChIPQC")
library(ChIPQC)
library(rtracklayer,quietly=TRUE)
#install.packages("DT")
library(DT)
library(dplyr)
library(tidyr)
#BiocManager::install("ChIPseeker")

library(ChIPseeker)
p<-arg_parser("Blacklist removal")
p<-add_argument(p, "blkList", help="black list bed file")
p<-add_argument(p, "peak_file",help="Open region peak file")
p<-add_argument(p, "bam_file", help="Open region bam file")
argv<-parse_args(p)

cat("QC for quality, duplicates and signal distribution")

blkList <- import.bed(argv$blkList)
openRegionPeaks <- argv$peak_file
qcRes <- ChIPQCsample(argv$bam_file, 
                      peaks = openRegionPeaks, annotation = "hg19", blacklist = blkList, 
                      verboseT = FALSE)

MacsCalls<-granges(qcRes)
data.frame(Blacklisted=sum(MacsCalls %over% blkList),
           Not_Blacklisted=sum(!MacsCalls %over% blkList))

MacsCalls_filtered<-MacsCalls[!MacsCalls %over% blkList]

