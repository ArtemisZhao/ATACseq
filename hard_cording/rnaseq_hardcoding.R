
library(Rsubread)

library(DESeq2)

info <- read.table('info.txt',header=TRUE)

####info$Bamfile leads to all the bam files
filenames <- file.path(paste0("/net/wonderland/home/zhaolinz/practice/RNAseq/bam/",info$Bamfile))

fc <- featureCounts(files=filenames,
                    
                    annot.inibuilt ="hg19",
                    
                    annot.ext= NULL,
                    
                    isPairedEnd=FALSE,
                    
                    nthreads=2)

colnames(fc$counts) <- row.names(coldata)

save(fc, file="countMatrix.RData")                       

cts <- (fc$counts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = info, design =  ~Psphen+Celltype+
                                Treat+skinhoming+Agebin+Celltype:Treat+
                                Celltype:skinhoming+skinhoming:Treat+Gender)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds <- DESeq(dds)

save(dds, file="RNA295_DESeq.RData")
