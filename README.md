# Running the ATAC-seq pipeline

## Fastq to BAM process

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


## Peak calling - Finding open regions

We use MACS2 in this step to call peak to the bam files. For detailed reference, please refer to https://github.com/taoliu/MACS.

```{r,eval=FALSE}
MACS2 callpeak -t ${bam_file}  -f BAMPE --outputdir $path_to_output\
--name ${output_name} -g hs --nomodel --shift -100 \
--extsize 200 -q 0.05\
```
By following this process, we would get 3 files.

- {output_name}_peaks.narrowPeak: Narrow peak format suitable for IGV and further analysis;

- {output_name}_peaks.xls: Peak table suitable for review in excel (not actually a xls but a tsv);

- {output_name}_summits.bed: Summit positions for peaks useful for finding motifs and plotting.


## Blacklisted regions removal

We apply bedtools to the original bam files to remove blacklisted regions. 

```{r,eval=FALSE}
bedtools intersect -v -abam ${bam_file} -b ${blklist.bed} > ${filtered_bam_file}
```

In this step, we will obtain the filtered bam files.

Notice that this serves as an alternative way of using ChIPQC to remove blacklisted regions. The guidance for the old one is as followed,

```{r,eval=FALSE}
Rscript blacklist_remove.R ${blkList.bed} ${peak_file} ${bam_file} 
```

## Differential expression ATAC-seq analysis

### DEseq2

#### 1. Identify non-redudant peaks

  - Define a set of non-redundant peaks present in at least 2 samples.

#### 2. Counting for differential ATAC-seq 

  - filter the peaks which are present in at least two replicates.

  - Use Rsubread to count paired reads landing in peaks.

#### 3. Contruct a DESeq2 object 

  - pass the count to __DESeqDataSetFromMatrix__ function so as to access these from DESeq2 later.

#### Command
```{r,eval=FALSE}
Rscript DEseq.R ${individualID} ${covariates} ${blkList.bed}\
--peak_dir {path_to_peakfiles} \
--bam_dir {path_to_bamfiles} \
--output_dir {path_to_output}
```

#### Input
While using the __DEseq.R__ script presented above, two text files that record the individual IDs and their corresponding conditions (covariates information) should be input, which refers to as {individualID} and {covariates} as bellow. 

- {covariates} should denote the treatment/control group that each bam file belongs to, e.g., 24h v.s. 0h.  

- {individualID} should record the ID of each file, and should be consistent with the sequence of the covariates. e.g., id1_24h, id2_0h... 

- {--peak_dir} points to a directory that preserves all the narrowPeak files which is output by the peak calling procedure. 

- Similarly, {--bam_dir} contains all the bam files that user intends to include in t he analysis, e.g., filtered bam files output by the blklist removal step.

#### Output

This step generates the following files for future usage, 

- {countMatrix}: count matrix records the counts landing in peaks;

- {atacDDS}: a DEseq object that is used for test any differences in ATAC-seq signal between groups.

#### Further analysis 

```{r}

```
### Limma 