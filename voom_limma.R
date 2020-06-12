library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mygene)
library(tracktables)
library(stringr)
library(limma)
library(edgeR)


load("0407.ATAC.countMatrix.RData")

coldata <- read.csv('0407.all.v2.csv',header=TRUE)

cts <- (fc$counts)

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





library(ChIPseeker)

library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(tracktables)

library(ReactomePA)

library(mygene)

library(clusterProfiler)



setwd('/net/wonderland/home/zhaolinz/practice/PCA')

#readPeak

DEG <- dir("/net/wonderland/home/zhaolinz/practice/PCA", pattern = "0407_anno.adj*", full.names = TRUE)

for(i in 1:length(DEG)){
  
  print(DEG[i])
  
  a <- read.table(DEG[i], header= F, fill= T)
  
  a = a[-1,]
  
  #dim(a) check if all the rows include
  
  #bed format covert to a Grange project
  
  a1 <- data.frame("chrom"=a$V2, "start"=a$V3, "end"= a$V4)
  
  df <-makeGRangesFromDataFrame(a1)
  
  #associatedGenes
  
  gene_DAR <- seq2gene(df, tssRegion = c(-3000, 3000), flankDistance = 100000, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, sameStrand = FALSE)
  
  geneName <- queryMany(gene_DAR, scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
  
  geneTable <- geneName$response
  
  gene <- data.frame("symbol" = geneTable$symbol,"EntrezID" = geneTable$query)
  
  ann.name <-substr(DEG[i],44,nchar(DEG[i])-0)
  
  genefile <-paste0("/net/wonderland/home/zhaolinz/practice/PCA/adjustSN0407/genetoDAR.",ann.name,sep="")
  
  write.table(gene, file=genefile, quote=F, sep="\t", row.names=F, col.names=T)
  
  #enrich.reactctome
  
  enrich.reactctome <- enrichPathway(gene_DAR, organism="human")
  
  enrich.reactctome<- as.data.frame(summary(enrich.reactctome ))
  
  reactomefile <- paste0("/net/wonderland/home/zhaolinz/practice/PCA/adjustSN0407/enrichRea.",ann.name,sep="")
  
  write.table(enrich.reactctome, file=reactomefile, quote=F, sep="\t", row.names=F, col.names=T)
  
  #enrich.GO.BP
  
  go = enrichGO(gene_DAR, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH", maxGSSize = 5000)
  
  go1 <-dropGO(go, level =1)
  
  go1 <-dropGO(go1, level =2)
  
  go1 <-dropGO(go1, level =3)
  
  go1 <-dropGO(go1, level =4)
  
  goBPfile <- paste0("/net/wonderland/home/zhaolinz/practice/PCA/adjustSN0407/enrichGO.5thBP.",ann.name,sep="")
  
  write.table(as.data.frame(go1),file=goBPfile, quote=FALSE, row.names=F,sep = "\t")
  
  #enrich.GO.MF
  
  go = enrichGO(gene_DAR, OrgDb = "org.Hs.eg.db", ont = "MF", pAdjustMethod = "BH", maxGSSize = 5000)
  
  go1 <-dropGO(go, level =1)
  
  go1 <-dropGO(go1, level =2)
  
  go1 <-dropGO(go1, level =3)
  
  go1 <-dropGO(go1, level =4)
  
  goMFfile <- paste0("/net/wonderland/home/zhaolinz/practice/PCA/adjustSN0407/enrichGO.5thMF.",ann.name,sep="")
  
  write.table(as.data.frame(go1),file=goMFfile, quote=FALSE, row.names=F,sep = "\t")
  
  #enrich.KEGG
  
  KEGG = enrichKEGG(gene_DAR, organism = "hsa", keyType = "kegg", pAdjustMethod = "BH")
  
  KEGGfile <-paste0("/net/wonderland/home/zhaolinz/practice/PCA/adjustSN0407/enrichKEGG.",ann.name,sep="")
  
  write.table(as.data.frame(KEGG),file=KEGGfile, quote=FALSE, row.names=F,sep = "\t")
  
}
