###load in differental peak region results
atacseqraw<-read.xlsx("~/Desktop/Pipeline/resultenrich/log2fccutSkinhoming_CLAp_vs_CLAn.xlsx")
atacseqraw<-atacseqraw[(abs(atacseqraw$log2FoldChange) >= 0.585),]

###load in DAR and DEG results
DAR<-read.table("~/Desktop/Pipeline/resultenrich/genetoDAR_0.585_Skinhoming_CLAp_vs_CLAn.txt",header=T)
DEG<-read.table("~/Desktop/Pipeline/rnaresult/RNAseq_skinhoming_CLAp_vs_CLAn.txt",header=T,fill=T)

DEG$geneid<-rownames(DEG)
DEG<-DEG[order(DEG$padj),]

#########1. contingency table summary
DEGlist<-unique(DEG$gene.symbol)
DARlist<-unique(DAR$symbol)
length(intersect(rownames(DEG),DAR$EntrezID))
length(union(rownames(DEG),DAR$EntrezID))


#########2. TSS distribution
maptoDEG<-atacseqraw[which(atacseqraw$gene.symbol %in% DEG$gene.symbol),]
maptonoDEG<-atacseqraw[-which(atacseqraw$gene.symbol %in% DEG$gene.symbol),]

quantile(maptoDEG$distanceToTSS,probs=c(0.05,0.95))
quantile(maptonoDEG$distanceToTSS,probs=c(0.05,0.95))

##histogram plot 
##please change the xlim accordingly
ggplot( data=maptoDEG,aes(x=distanceToTSS)) +
  geom_histogram( bins=100, alpha=0.9) +
  xlim(-147178,143768)+
  ylim(0,1500)+
  ggtitle("TSS histogram for chrom region map to DEGs") +
  theme(plot.title = element_text(size=15))


ggplot( data=maptonoDEG,aes(x=distanceToTSS)) +
  geom_histogram(bins =100 , alpha=0.9) +
  xlim(-246869,338149)+
  ggtitle("TSS histogram for chrom region map to non DEGs") +
  theme(plot.title = element_text(size=15))

##density plot
ggplot( data=data.frame(value=c(maptoDEG$distanceToTSS,maptonoDEG$distanceToTSS),
                        label=c(rep("map to DEG",nrow(maptoDEG)),rep("map to nonDEG",nrow(maptonoDEG)))
                        ),aes(x=value,color=label)) +
  geom_density(alpha=0.9) +
  xlim(-246869,338149)+
  ggtitle("Density plot for TSS distance") +
  xlab("TSS distance")+
  theme(plot.title = element_text(size=15))+
  scale_color_brewer(palette="Dark2") + theme_minimal()
ggplot( data=data.frame(value=c(maptoDEG$distanceToTSS,maptonoDEG$distanceToTSS),
                        label=c(rep("map to DEG",nrow(maptoDEG)),rep("map to nonDEG",nrow(maptonoDEG)))
),aes(x=label,y=value,color=label)) +
  geom_boxplot()+
  ggtitle("Box plot for TSS distance") +
  ylab("TSS distance")+
  ylim(-147178,143768)+
  theme(plot.title = element_text(size=15))+
  scale_color_brewer(palette="Dark2") + theme_minimal()


#########3. rank correlation 
DEGinter<-DEG[which(rownames(DEG) %in% intersect(rownames(DEG),DAR$ID)),]
DARinter<-DAR[which(DAR$ID %in% intersect(rownames(DEG),DAR$ID)),]

DEGinter$rank<-c(1:1678)
DEGinter$ID<-rownames(DEGinter)
DARinter$rank<-c(1:1678)

rankplot<-left_join(DEGinter,DARinter,by="ID")


