library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mygene)
library(tracktables)
library(stringr)
library(limma)
library(edgeR)

#setwd('/net/wonderland/home/zhaolinz/practice/PCA')

#load("0407.ATAC.countMatrix.RData")
p<-arg_parser("Differential expression analysis")
p<-add_argument(p,"count_matrix", help="Matrix that records count landing in peaks")
p<-add_argument(p,"",help=)
argv<-parse_args(p)

coldata <- read.csv('0407.all.v2.csv',header=TRUE)

cm<-argv$count_matrix

cts <- (cm$counts)

colnames(cts) <- coldata$SampleID

dim(cts)

d0 <- DGEList(cts)

d0 <- calcNormFactors(d0) #TMM default,the weighted trimmed mean of M-values (to the reference)

#to remove the low reads region

max(rowSums(d0$counts))

min(rowSums(d0$counts))

dim(d0[rowSums(d0$counts)> 6960,]) #average reads per lib >10

d <- d0[rowSums(d0$counts)> 6960,]

dim(d)

f <- factor(coldata$Celltype.Treats)

# cell means model with Celltype.Treats, AverSN as blocking factors

design <- model.matrix(~0 + f + AverSN, coldata)

head(design)

y <- voomWithQualityWeights(d, design= design, plot = T)

dev.off()

save(y, file="0407.ATAC.voomWithQualityWeights.RData")

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

fit3  <- eBayes(fit2, trend = TRUE)



#plotSA(fit)
#Plot residual standard deviation versus average log expression for a fitted microarray linear model.

#dev.off()

DAR.ID <- c("CD4CLA.0h","CD8CLA.0h", "CD4CLA.24h", "CD8CLA.24h", 
            "CD4CLAN.ST","CD4CLAP.ST", "CD8CLAN.ST","CD8CLAP.ST")

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



