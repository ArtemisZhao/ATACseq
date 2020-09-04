
library(Rsubread)

library(DESeq2)

library(ChIPseeker)

library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(tracktables)

library(ReactomePA)

library(mygene)

library(clusterProfiler)

library(DOSE)

info <- read.table('rnainfo.txt',header=TRUE)

####info$Bamfile leads to all the bam files
filenames <- file.path(paste0("/net/wonderland/home/zhaolinz/practice/RNAseq/bam/",info$bamlist))

fc <- featureCounts(files=filenames,
                    
                    annot.inbuilt ="hg19",
                    
                    annot.ext= NULL,
                    
                    isPairedEnd=FALSE,
                    
                    nthreads=2)

colnames(fc$counts) <- row.names(info)

save(fc, file="countMatrix.RData")                       

cts <- (fc$counts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = info, design =  ~Psphen+Celltype+
                                Treat+skinhoming+Celltype:Treat+
                                Celltype:skinhoming+skinhoming:Treat+Gender)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)

save(dds, file="RNA295_DESeq.RData")

###LRT for skinhoming
#reduce_design<-as.formula(~Psphen+Celltype+
                           # Treat+Agebin+Celltype:Treat+Gender)

#dds_lrt <- DESeq(dds, test="LRT", reduced = reduce_design)

#lrt_skinhoming<-results(dds_lrt)
#write.table(lrt_skinhoming, file="./lrt_skinhoming.txt", quote=F, sep="\t", row.names=T, col.names=T)


####write results
resultname1<-resultsNames(dds)

###Or extract what you want
resultname <- resultname1[-1]


###DEG and write results

a <- lapply(resultname,
            function(x){
              res <- results(dds, name=x) 
              res <- res[order(res$pvalue),]
              res <- res[(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]
              geneName<-queryMany(rownames(res),scopes="entrezgene", fields="symbol", species="human",returnall=TRUE)
              res$gene.symbol<-geneName$response$symbol
              res
            })

for (i in 1:length(a)){
  write.table(a[[i]], file=paste0("./",resultname1[i],".txt"), quote=F, sep="\t", row.names=T, col.names=T)
}



#####enrichGO and enrichKEGG by OSA
for (i in 1:length(a)){
  
  enrich.reactctome <- enrichPathway(rownames(a[[i]]), organism="human",readable=TRUE)
  
  enrich.reactctome<- as.data.frame(summary(enrich.reactctome ))
  
  write.table(enrich.reactctome, file=paste0("./enrichreact_",resultname[i],".txt"), quote=F, sep="\t", row.names=F, col.names=T)
  
  #enrich.GO.BP
  
  go = enrichGO(rownames(a[[i]]), OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH", maxGSSize = 5000,readable=TRUE)
  
  go1 <-dropGO(go, level =1)
  
  go1 <-dropGO(go1, level =2)
  
  go1 <-dropGO(go1, level =3)
  
  go1 <-dropGO(go1, level =4)
  
  write.table(as.data.frame(go1),file=paste0("./enrichGO_BP_",resultname[i],".txt"), quote=F, row.names=F,sep = "\t")
  
  #enrich.GO.MF
  
  go = enrichGO(rownames(a[[i]]), OrgDb = "org.Hs.eg.db", ont = "MF", pAdjustMethod = "BH", maxGSSize = 5000,readable=TRUE)
  
  go1 <-dropGO(go, level =1)
  
  go1 <-dropGO(go1, level =2)
  
  go1 <-dropGO(go1, level =3)
  
  go1 <-dropGO(go1, level =4)
  
  write.table(as.data.frame(go1),file=paste0("./enrichGO_MF_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
  
  #enrich.KEGG
  
  KEGG = enrichKEGG(rownames(a[[i]]), organism = "hsa", keyType = "kegg", pAdjustMethod = "BH")
  
  y <- setReadable(KEGG, 'org.Hs.eg.db', keyType="ENTREZID")
  
  write.table(as.data.frame(y),file=paste0("./enrichKEGG_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
}


##### enrich analysis by GSEA

b <- lapply(resultname,
            function(x){
              res <- results(dds, name=x) 
              res
            })

for(i in 1:length(b)){
  
  midterm<-b[[i]]
  genelist=midterm[ ,2]
  
  names(genelist)= rownames(midterm)
  
  genelist= sort(genelist,decreasing=T)
  
  gse.reactome <-gsePathway(genelist, organism="human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  gse.reactome<-setReadable(gse.reactome,'org.Hs.eg.db', keyType="ENTREZID")
  
  write.table(summary(gse.reactome), file=paste0("./gsereactome_",resultname[i],".txt"), quote=F, sep="\t", row.names=F, col.names=T)
  
  go=gseGO(genelist, OrgDb= "org.Hs.eg.db", ont= "BP", nPerm= 100000,minGSSize= 100, maxGSSize=500, pvalueCutoff =0.05, pAdjustMethod = "BH")
  
  go <-setReadable(go,'org.Hs.eg.db', keyType="ENTREZID")
  
  write.table(summary(go),file=paste0("./gsego_bp_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
  
  go1=gseGO(genelist, OrgDb= "org.Hs.eg.db", ont= "MF", nPerm= 100000,minGSSize= 100, maxGSSize=500, pvalueCutoff =0.05, pAdjustMethod = "BH")
  
  go1<-setReadable(go1,'org.Hs.eg.db', keyType="ENTREZID")
  
  write.table(summary(go1),file=paste0("./gsego_mf_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
  
  KEGG = gseKEGG(genelist, organism = "hsa", pvalueCutoff = 0.05, keyType = "kegg", pAdjustMethod = "BH")
  
  KEGG <-setReadable(KEGG,'org.Hs.eg.db', keyType="ENTREZID")
  
  write.table(summary(KEGG),file=paste0("./gsekegg_",resultname[i],".txt"), quote=FALSE, row.names=F,sep = "\t")
  
}
