###YIZHAO ATAC-seq Fastq to BAM processing
#install.packages("BiocManager")
library(BiocManager,quietly = TRUE)
#BiocManager::install("Rsubread")
library(Rsubread,quietly = TRUE)
#BiocManager::install("Rsamtools")
library(Rsamtools,quietly = TRUE)
#BiocManager::install("GenomicAlignments")
library(GenomicAlignments,quietly = TRUE)
#BiocManager::install("rtracklayer")
library(rtracklayer,quietly = TRUE)
#install.packages("ggplot2")
library(ggplot2,quietly = TRUE)
#install.packages("magrittr")
library(magrittr,quietly = TRUE)
#install.packages("argparser")
library(argparser,quietly = TRUE)
#BiocManager::install("GenomicAlignments")
library(GenomicAlignments,quietly=TRUE)

p<-arg_parser("ATAC-seq fastq to BAM processing")

p<-add_argument(p, "genome_file", help="Reference genome fasta file")
p<-add_argument(p, "input_file_1", help="Sequence reads file to be aligned; first reads file for paired-end reads")
p<-add_argument(p, "input_file_2",  help="Second reads file for paired-end reads")
p<-add_argument(p, "prefix", help="Prefix for output bam files; usually <sample_id>")
p<-add_argument(p, "--seqnames",default=NULL, help="Reads mapping to specific seqnames")
p<-add_argument(p, "--output_dir", short="-o", default=".", help='Output directory')
p<-add_argument(p, "--maxFragLength",default=2000, help="Maximum fragment length in alignment")
argv<-parse_args(p)

cat("Generating a reference genome ...")

indexForSubread <- paste0(dirname(argv$genome_file),gsub("\\.fa$", "", argv$genome_file))

buildindex(indexForSubread, argv$genome_file, indexSplit = FALSE)

#############################################################
cat("Aligning reads to reference genome ...")

outBAM<-paste0(argv$output_dir,"/",argv$prefix,".bam",sep="")

align(indexForSubread, readfile1 = argv$input_file_1, readfile2 = argv$input_file_2, 
      output_file = outBAM, nthreads = 2, type = 1, unique = TRUE, 
      maxFragLength = as.numeric(argv$maxFragLength))

#############################################################
cat("Sorting and indexing ...")

sortedBAM<-file.path(dirname(outBAM),paste0("Sorted_",argv$prefix,".bam",sep=""))

sortBam(outBAM, gsub("\\.bam","",basename(sortedBAM)))
indexBam(sortedBAM)

###############################################################
cat("Calculating mapped reads ...")
pmapped<-propmapped(sortedBAM)
print(pmapped)

################################################################

if (!is.null(argv$seqnames)){
cat("Post-alignment processing ...")
cat(paste0("Reading in mapped reads to",argv$seqnames))

atacReads<-readGAlignmentPairs(sortedBAM,param=ScanBamParam(mapqFilter = 1,
          flag=scanBamFlag(isPaired = TRUE, isProperPair = TRUE),what=c("qname",
          "mapq","isize"), which =GRanges(argv$seqnames,IRanges(1,63025520))))

##################################################################
cat("Creating BAM files split by insert sizes. ...")
atacReads_Open <- atacReads[insertSizes < 100, ]
atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]

openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM)
monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM)
diNucBam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM)

export(atacReads_Open, openRegionBam, format = "bam")
export(atacReads_MonoNuc, monoNucBam, format = "bam")
export(atacReads_diNuc,diNucBam,format = 'bam')
}
