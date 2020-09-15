library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mygene)
library(tracktables)
library(stringr)
library(limma)
library(edgeR)



p<-arg_parser("Differential expression analysis - limma combined with voom")
p<-add_argument(p,"peak_type", help="Table contains peak type names for each file")
p<-add_argument(p, "group_name",help="Table contains covariate factor")
p<-add_argument(p, "blkList", help="black list bed file")
p<-add_argument(p, "--peak_dir", short="-pdir", default=".", help='Peak file directory')
p<-add_argument(p,"--bam_dir",short="-bdir", default=".", help='Bam file directory')
p<-add_argument(p, "--output_dir",short="-o",default=".",help="Output directory")
p<-add_argument(p, "--add_cov", default=NULL, help="Additional covariate, e.g. signal to noise ratio")
p<-add_argument(p,"--cont",help="Whether the additional covariate is continuous",flag=TRUE)
p<-add_argument(p,"--filter",help="The lower bound of the counts to be filtered in voom",default=6900)
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
  if (isTRUE(argv$cont)){
    AddGroup <- as.numeric(t(read.table(argv$add_cov)))
    coldata<-data.frame(Group=Group,AddGroup=AddGroup)
  }else{
    AddGroup<-factor(t(read.table(argv$add_cov)))
    coldata<-data.frame(Group=Group,AddGroup=AddGroup)
    }
}else{coldata<-data.frame(Group=Group)}

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

cm<-myCounts

cts <- (cm$counts)

dim(cts)

d0 <- DGEList(cts)

d0 <- calcNormFactors(d0) 

max(rowSums(d0$counts))

min(rowSums(d0$counts))

cat("Filtering the low counts data...\n")
dim(d0[rowSums(d0$counts)> as.numeric(argv$filter),]) #average reads per lib >10

d <- d0[rowSums(d0$counts)> as.numeric(argv$filter),]

dim(d)

# cell means model with Celltype.Treats, AverSN as blocking factors


if (!argv$add_cov=="NULL"){
  design <- model.matrix(~0 + Group + AddGroup, coldata)
}else{design<-model.matrix(~0+Group,coldata)}

head(design)

cat("Constructing voom weights...\n")
pdf(file=paste0(argv$output_dir,"voomweight.pdf"),width=11,height=6,bg="white",pointsize=14)
y <- voomWithQualityWeights(d, design= design, plot = T)
dev.off()

save(y, file="voomWithQualityWeights.RData")

head(cpm(y, log=F))#to yield the normalized counts from y.

fit <- lmFit(y, design)


cont.matrix <- makeContrasts(P1="f0hctCD4CLAP-f0hctCD4CLAN",
                             P2="f0hctCD8CLAP-f0hctCD8CLAN", 
                             P3="f24hstCD4CLAP-f24hstCD4CLAN", 
                             P4="f24hstCD8CLAP-f24hstCD8CLAN",
                             P5="f24hstCD4CLAN-f0hctCD4CLAN", 
                             P6="f24hstCD4CLAP-f0hctCD4CLAP",
                             P7="f24hstCD8CLAN-f0hctCD8CLAN",
                             P8="f24hstCD8CLAP-f0hctCD8CLAP", levels = design)

fit2  <- contrasts.fit(fit, cont.matrix)

cat("Fitting limma model...\n")
fit3  <- eBayes(fit2, trend = TRUE)


DAR.ID <- c("CD4CLA.0h","CD8CLA.0h", "CD4CLA.24h", "CD8CLA.24h", 
            "CD4CLAN.ST","CD4CLAP.ST", "CD8CLAN.ST","CD8CLAP.ST")

cat("Saving DARs...\n")
for(i in 1:8){
  
  res <- topTable(fit3,coef=i,n = Inf,adjust.method="BH",sort.by="p")
  
  res <- res[(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & abs(res$logFC) >= 0.585), ]
  
  length(which(res$adj.P.Val < 0.05))
  
  #convert dataframe to grange project
  
  res_split <- str_split_fixed(rownames(res), "_", 4)
  
  res$Chrom <- res_split[ ,2]
  
  res$Start <- res_split[ ,3]
  
  res$End <- res_split[ ,4]
  
  res_grange <- makeGRangesFromDataFrame(res,keep.extra.columns=TRUE,seqnames.field="Chrom",start.field="Start",end.field="End",starts.in.df.are.0based=TRUE)
  
  anno_res <- annotatePeak(res_grange, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  anno_res <-as.GRanges(anno_res)
  
  gr <- as(anno_res , "data.frame")
  
  GeneName <- queryMany(gr$geneId, scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
  
  gr$gene.symbol <- GeneName$response$symbol
  
  file.ID <- paste0("0407_anno.adjSN.",DAR.ID[i],".DAR.txt",sep="")
  
  write.table(gr,file=file.ID, quote=F, sep="\t", row.names=T, col.names=T)
  
}



