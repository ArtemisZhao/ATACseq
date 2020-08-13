########### DESeq2 Likelihood Ratio Test --- HARD-CODING
########### YI ZHAO
library(ChIPQC)
library(ChIPseeker)
library(DT)
library(tidyr)
library(Rsubread)
library(dplyr)
library(rtracklayer)
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

####load countmatrix
myCounts<-load("./countMatrix.RData")

####load ConsensusToCount
consensusToCount<-load("./consensusToCount.Rdata")

####load covariate information
fullinfo<-read.table("./fullinfo.txt")

####check the correctness of the variable names
names(fullinfo)

###full model: please make sure the model specification is correct
atacDDS <- DESeqDataSetFromMatrix(myCounts, fullinfo, ~Psphen+Age+Gender+SN+Time+Skinhoming+Celltype+
                                    Time:Skinhoming+Time:Celltype+Skinhoming:Celltype, 
                                  rowRanges = consensusToCount) 

###reduced model : for Skinhoming variable
reduce_design<-as.formula(~Psphen+AgeGender+SN+Time+Celltype+Time:Celltype)

dds_lrt <- DESeq(atacDDS, test="LRT", reduced = reduce_design)

lrt_skinhoming<-results(dds_lrt)


###reduced model : for stimulation variable
reduce_design<-as.formula(~Psphen+AgeGender+SN+Skinhoming+Celltype+Skinhoming:Celltype)

dds_lrt <- DESeq(atacDDS, test="LRT", reduced = reduce_design)

lrt_stimulation<-results(dds_lrt)

###write results
write.table(lrt_skinhoming, file="./lrt_skinhoming.txt", quote=F, sep="\t", row.names=T, col.names=T)

#write.table(lrt_stimulation, file="./lrt_stimulation.txt", quote=F, sep="\t", row.names=T, col.names=T)
