# ATACseq
# Requirements
We will use Bioconductor package in R to accomplish the following ATACseq steps. For more installation details, please refer to https://www.bioconductor.org/install/.

The updated version of Bioconductor is 3.10, which requires R with the version >=3.6.

```{r,eval=FALSE}
install.packages("BiocManager")
library(BiocManager)
```

Using __BiocManager::install()__ is the recommended way to install Bioconductor packages; we will show how and what specific packages to installed in the relevent steps.

# Fastq to BAM process
We will process the dataset from fastq format to BAM for further analysis.

## Generate the reference genome

```{r,eval=FALSE}
BioManager::install(Rsubread)
library(Rsubread)

genome <- "$PathToFasea/Homo_sapiens_assembly19.fasta "
indexForSubread <- gsub("\\.fa$", "", genome)

buildindex(indexForSubread, genome, indexSplit = FALSE)
#save(indexForSubread,file="indexForSubread.Rdata")
```

The reference genome with fasta format is required in the step, which is located under the $PathToFasta directory.

The generated reference is stored in the R subject named as __indexForSubread__.

If the user will perform the following aligning step immediately, the last step of storing index could be skipped.

## Align reads to reference genome

In this step, we align the paired-end ATAC-seqp reads using the Rsubread function.

```{r,eval=FALSE}
library(Rsubread)
read1<-"$PathToData/#read1name.fastq.gz"
read2<-"$PathToData/#read2name.fastq.gz"
outBAM<-"$PathtoSave/#outputname.bam"

align(indexForSubread, readfile1 = read1, readfile2 = read2, output_file = outBAM, 
    nthreads = 2, type = 1, unique = TRUE, maxFragLength = #maxvalue))
```

In this step, we require the user to provide two reads data in __fastq.gz__ format for excuting aligning. Notice that they are stored under the $PathToData directory and the name #read1name and #read2name should be modified to their name correspondingly.

The name for the output bam file #outputname is user-defined. The bam file will be instored under $PathToSave directory.

Parameter __maxFragLength__ is also a user-defined value, it should be set based on the study which generates the data file.

## Sort and index
In this step, we sort and index the output BAM file in the previous step for use in external tool.

First we sort our aligned data by sequence order.

We then index the sorted bam file by Rsamtools package available in Bioconductor.

```{r,eval=FALSE}
BiocManager::install("Rsamtools")
library(Rsamtools)

#sort
sortedBAM <- file.path(dirname(outBAM), paste0("Sorted_", basename(outBAM)))
sortBam(outBAM, gsub("\\.bam", "", basename(sortedBAM)))

#index
indexBam(sortedBAM)
```

The sorted and indexed bam file is stored under the same directory with the output bam file. And the name for sorted bam file is __Sorted_#outputname.bam__. 

In order to get access to the sorted bam file generated in this step, we use the following command. (Replace the $PathToSave and #outputname with the ones defined in the aligning step.)

```{r,eval=FALSE}
sortedBAM<-"$PathToSave/Sorted_#outputname.bam"
```

## Summarize and visualize 
This step provides the summary(count) and the distribution plot for the mapped reads, which is optional for excution.

- Number of mapped reads
```{r}
library(Rsubread)
#sortedBAM<-"$PathToSave/Sorted_#outputname.bam"
nummapped <- propmapped(sortedBAM)
print(nummapped)
```

- Distribution of mapped reads
```{r}
library(Rsamtools)
require(ggplot2)
require(magrittr)

idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
```


# Post-alignment process
In the following stwps, we will work with alignments.

## Read in mapped reads

We read our aligned data by using the GenomicAlignments package and save in the __atacReads__ element.

```{r}
BiocManager::install("GenomicAlignments")
library(GenomicAlignments)

#sortedBAM<-"$PathToSave/Sorted_#outputname.bam"
atacReads <- readGAlignmentPairs(sortedBAM, param = ScanBamParam(mapqFilter = 1, 
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
        "mapq", "isize"), which = GRanges("chr20", IRanges(1, 63025520))))

atacReads
#save(atacReads, file="atacReads.Rdata")
```

We read properly paired, uniquely mapped reads in this step. Acting as an example, we pull reads that map to Chromosome 20 (which is specifed in __Granges("chr20",...)__, users can define their own interest). 

We also retrieve the read name, mapq scores and the insert sizes, specified as __c("qname", "mapq", "isize")__.

## Retrieve insert sizes
The paired aligned data is read into R element __acatReads__ we can retreive the insert sizes from the elementMetadata attached to reads. 

```{r}
library(GenomicAlignments)
#load("atacReads.Rdata")

atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)
```

We show the head of insertSized for read1 in the step.

## Subsetting reads by insert sizes and create coresponding bam files
In this step, we divide the aligned reads into reads representing nucleosome free and nucleosome occupied.

Here we create BAM files for the reads representing nucleosome free, mono and di nucleosome by using insert sizes to filtering read and write them back to BAM files 

```{r}
atacReads_Open <- atacReads[insertSizes < 100, ]
atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]

openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM)
monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM)
diNucBam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM)

BiocManager::install("rtracklayer")
library(rtracklayer)
export(atacReads_Open, openRegionBam, format = "bam")
export(atacReads_MonoNuc, monoNucBam, format = "bam")
export(atacReads_diNuc,diNucBam,format = 'bam')
```

Notice that the files will be named after __Sorted_#outputname_openRegions.bam__, __Sorted_#outputname_monoNuc.bam__ and __Sorted_#outputname_diNuc.bam__, and exported to the current directory.

# Find open regions
## Peak calling using MACS2
As we discussed, the reads with insert sizes less than 100bp represent fragments coming from open chromatin and around transcription factors bound to the DNA.

There are many methods available to call nucleosome free regions from ATAC-seq data with many borrowed from ChIP-seq analysis.

We will use __MACS2__ in the following steps. For detailed reference, please see https://github.com/taoliu/MACS.

##Single end peak calling
With single end sequencing from ATAC-seq we do not know how long the fragments are.

To identify open regions therefore requires some different parameters for MACS2 peak calling.

One strategy employed is to shift read 5’ ends by -100 and then extend from this by 200bp. Considering the expected size of our nucleosome free fragments this should provide a pile-up over nucelosome regions suitable for MACS2 window size.
```{r,eval=FALSE}
MACS2 callpeak -t singleEnd.bam --nomodel --shift -100 --extsize 200 
--format BAM -g MyGenome
```

Alternatively for the nucleosome occupied data we can adjust shift and extension to centre the signal on nucleosome centres (nucleosomes wrapped in 147bp of DNA).

```{r,eval=FALSE}
MACS2 callpeak -t singleEnd.bam --nomodel --shift 37 --extsize 73 
--format BAM -g MyGenome
```

## Paired end peak calling
If we have sequenced paired-end data then we do know the fragment lengths and can provide pre-filter BAM files to MACS2.

We have to simply tell MACS2 that the data is paired using the format argument. By default MACS2 will guess it is single end BAM.

(Apply on sorted bam is not proper and not what the user is interested, just for a reference.)
```{r, eval=FALSE}
$PathToMACS2/macs2 callpeak 
-t $PathToSave/Sorted_#outputname.bam 
-f BAMPE 
--outdir $PathToSave 
--name #NAME -g hs
```

Following peak calling we would get 3 files.

- NAME.narrowPeak: Narrow peak format suitable for IGV and further analysis

- NAME_peaks.xls: Peak table suitable for review in excel.(not actually a xls but a tsv)

- summits.bed: Summit positions for peaks useful for finding motifs and plotting

## Blacklisted regions removal --- ChIPQC

```{r,eval=FALSE}
Rscript blacklist_remove.R ${blkList.bed} ${peak_file} ${bam_file} 
```
