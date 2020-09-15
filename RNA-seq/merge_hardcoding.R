remainbam<-read.xls("Yi_482library information.xlsx",header=TRUE,sheet=1)

a<-c()
for (i in 1:187){
  temp <-gsub("_","",remainrnabam$X.RNAseq_RunID.[i])
  runid<-unlist(strsplit(temp,", "))
  
  coreid<-unlist(strsplit(as.character(remainrnabam$RNAseq_CoreID[i]),", "))
  
  ### if you have 4 bam files please add in one loop
  if (length(runid)==3){
    a[i]<-paste0("samtools merge ",remainrnabam$bamlist[i],
                 " /net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[1],
                 "/Sample_",coreid[1],".merged.bam ",
                 "/net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[2],
                 "/Sample_",coreid[2],".merged.bam ",
                 "/net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[3],
                 "/Sample_",coreid[3],".merged.bam",sep="")
  }
  
  if (length(runid)==2){
    a[i]<-paste0("samtools merge ",remainrnabam$bamlist[i],
                 " /net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[1],
                 "/Sample_",coreid[1],".merged.bam ",
                 "/net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[2],
                 "/Sample_",coreid[2],".merged.bam",sep="")
  }
  
  if (length(runid)==1){
    a[i]<-paste0("samtools merge ",remainrnabam$bamlist[i],
                 " /net/psorgen/alextsoi/Psoriasis_Tcells_RNAseq/STAR/*",runid[1],
                 "/Sample_",coreid[1],".merged.bam ")
  }
  
}

write.table(a,file="command.txt",row.names = F,col.names = F,quote=F)


