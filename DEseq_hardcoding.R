########### DESeq2 --- HARD-CODING VERSION
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

####hard-cording part, please make sure the model specification is correct
atacDDS <- DESeqDataSetFromMatrix(myCounts, fullinfo, ~Psphen+Age+Gender+SN+Time+Skinhoming+Celltype+
                                    Time:Skinhoming+Time:Celltype+Skinhoming:Celltype, 
                                  rowRanges = consensusToCount) 
atacDDS <- DESeq(atacDDS)

resultname<-resultsNames(atacDDS)

####Extract the result you want: 
####Time, Skinhoming, Celltype, Time:Skinhoming, Time:Celltype, Skinhoming:Celltype
resultname<-resultname[c(6,7,8,9,10,11)]

a <- lapply(resultname,
            function(x){
              print(x)
              res <- results(atacDDS, name=x,format = "GRanges") 
              res <- res[order(res$pvalue)]
              res <- res[(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]
              gr <- annotatePeak(res,TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
              gr<-as.data.frame(as.GRanges(gr))
              GeneName <- queryMany(gr$geneId, scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
              gr$gene.symbol <- GeneName$response$symbol
              gr
            })

####write the outcome to table
for (i in 1:length(a)){
  write.table(a[[i]], file=paste0("./",resultname[i],".txt"), quote=F, sep="\t", row.names=T, col.names=T)
}


####Can also extract information as the following
#### Stimulation effect in CD8 cell
#b<-results(dds, list( c("Time_24h_vs_0h","Time24h.CelltypeCD8") ),format = "GRanges") 


### Enrichment analysis

for (i in 1:length(a)){
 geneName <- queryMany(a[[i]], scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
 geneTable <- geneName$response
 gene <- data.frame("symbol" = geneTable$symbol,"EntrezID" = geneTable$query)
 write.table(gene, file=paste0("./genetoDAR_",resultname[i],".txt"), quote=F, sep="\t", row.names=F, col.names=T)

 #entich.reacrctome
 enrich.reactctome <- enrichPathway(a[[i]], organism="human",readable=TRUE)
 enrich.reactctome<- as.data.frame(summary(enrich.reactctome ))
 write.table(enrich.reactctome, file=paste0("./reactome_",resultname[i],".txt"), quote=F, sep="\t", row.names=F, col.names=T)


 #enrich.GO.BP
 go = enrichGO(a[[i]], OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH", maxGSSize = 5000,readable=TRUE)
 go1 <-dropGO(go, level =1)
 go1 <-dropGO(go1, level =2)
 go1 <-dropGO(go1, level =3)
 go1 <-dropGO(go1, level =4)
 write.table(as.data.frame(go1),file=paste0("./enrichGO.5thBP_",resultname[i],".txt"), quote=F, row.names=F,sep = "\t")

#enrich.GO.MF
 go = enrichGO(a[[i]], OrgDb = "org.Hs.eg.db", ont = "MF", pAdjustMethod = "BH", maxGSSize = 5000,readable=TRUE)
 go1 <-dropGO(go, level =1)
 go1 <-dropGO(go1, level =2)
 go1 <-dropGO(go1, level =3)
 go1 <-dropGO(go1, level =4)
 write.table(as.data.frame(go1),file=paste0("./enrichGO.5thMF_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")

#enrich.KEGG
 KEGG = enrichKEGG(a[[i]], organism = "hsa", keyType = "kegg", pAdjustMethod = "BH")

 y <- setReadable(KEGG, 'org.Hs.eg.db', keyType="ENTREZID")

 write.table(as.data.frame(y),file=paste0("./enrichKEGG_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
}




