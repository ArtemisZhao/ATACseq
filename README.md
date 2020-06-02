## Running the pipeline

### Fastq to BAM process

1. Generate the reference genome

2. Align reads to the reference genome

3. Sort and index

4. Read in mapped reads

5. Subsetting reads by insert sizes

6. Visualize(*)

```{r,eval=FALSE}
Rscript fastqtobam.R {genome_file} {input_1} {input_2} {prefix}\
 --seqnames chr20\
 --output_dir ./output\
 --maxFragLength 2000 
```

### Peak calling - Finding open regions

We will use MACS2 in the following steps. For detailed reference, please see https://github.com/taoliu/MACS.

```{r,eval=FALSE}
MACS2 callpeak 
-t $outputbam.bam\
-f BAMPE \
--outputdir ./output\
--name $NAME \
-g hs
```

## Blacklist removal

1. QC for low quality, duplicates and signal distribution

2. Remove blacklisted peaks 

```{r, eval=FALSE}
Rscript blacklist_remove.R {blkList.bed} {peak_file} {bam_file}
```

## Differential expression ATAC-seq analysis

1. Identify non-redudant peaks

2. Counting for differential ATAC-seq

3. DESeq2 for differential ATAC-seq

```{r,eval=FALSE}
Rscript DESeq.R {peak_type} {group} {blkList.bed}\
--peak_dir ./peakfiles \
--bam_dir ./bamfiles \
--output_dir ./output
```
