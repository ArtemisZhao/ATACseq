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

This process generates 3 files,

- {output_name}_peaks.narrowPeak: Narrow peak format suitable for IGV and further analysis;

- {output_name}_peaks.xls: Peak table suitable for review in excel (not actually a xls but a tsv);

- {output_name}_summits.bed: Summit positions for peaks useful for finding motifs and plotting.


## Blacklisted regions removal

We apply bedtools to the original bam files to remove blacklisted regions. 

```{r,eval=FALSE}
bedtools intersect -v -abam ${bam_file} -b ${blklist.bed} > ${filtered_bam_file}
```

In this step, we will obtain the filtered bam files named as {filtered_bam_file}.

Notice that this serves as an alternative way of using ChIPQC to remove blacklisted regions. The guidance for the old one can be found in the oldreadme.


## Differential expression ATAC-seq analysis

### I. DEseq2

#### Steps 

 1. Identify non-redudant peaks 
 

 2. Counting for differential ATAC-seq - use __Rsubread__ to count paired reads landing in peaks.

 3. Contruct a DESeq2 object - Pass the count to __DESeqDataSetFromMatrix__ function 
 
 
 4. Annotate peaks to genes
 

#### Command
```{r,eval=FALSE}
Rscript DEseq.R ${ID} ${covariate} ${blkList.bed}\
--peak_dir {path_to_peakfiles} \
--bam_dir {path_to_bamfiles} \
--output_dir {path_to_output} \
```

#### Input
While using the __DEseq.R__ script presented above, two text files that record the IDs and the corresponding conditions (covariates information) should be input, which refers to as {individualID} and {covariates} as above. 


- {ID} should record the ID of each file, and should be consistent with the sequence of the covariate information. e.g., id1_24h, id2_0h... 

- {covariate} contains the covariate information;

- {--peak_dir} points to a directory that preserves all the narrowPeak files output by the Peak calling step. 

- {--bam_dir} contains all the bam files that user intends to include in the analysis, e.g., filtered bam files output by the blacklist removal step.



#### Output

The DEseq2 step generates the following files for future usage and summarize, 

- {countMatrix}: count matrix records the counts landing in peaks;

- {atacDDS}: a DEseq object that is used for test any differences in ATAC-seq signal between groups.

- Files that records the contrast(group) information, the corresponding p-values for the siginificant peaks, and genes anotated to the peaks. These files are named by following the pattern of {Group_CD4CLAp0hct_vs_CD4CLAn0hct.txt}, in which CD4CLAp0hct anc CD4CLAn0hct are replaced by different groups existing in the covariate factor.


### II. Limma and voom

#### Steps 

 1. Identify non-redudant peaks 

 2. Counting for differential ATAC-seq - use __Rsubread__ to count paired reads landing in peaks.

 3. Construct voom weight.
 
 4. Perform limma model on different contrasts.
 
 5. Annotate peaks to nearest genes.
 
 6. Save DAR files to output directory

#### Command
```{r,eval=FALSE}
Rscript voom_limma.R ${ID} ${covariate} ${blkList.bed}\
--peak_dir {path_to_peakfiles} \
--bam_dir {path_to_bamfiles} \
--output_dir {path_to_output} \
--add_cov ${StN_ratio}\
--filter 6900
```

#### Input
While using the __limma_voom.R__ script presented above, please refer to __DEseq.R__ guidance, they follow the same input rules.


- {--filter} determines the lower bound of counts to be included while calculating the voom weight. Notice the output __voomweight.pdf__ plot should be utilized to determine whether the {--filter} parameter is well defined.

#### Output

The voom_limma step generates the following files for future usage and summarize, 

- {countMatrix}: count matrix records the counts landing in peaks;

- {voomweight.pdf} : plot presented voom weight ( mean versus variance relationship; and weight for each individual)

- {voomWithQualityWeights.RData}

- Files that records the contrast(group) information, the corresponding p-values for the siginificant peaks, and genes anotated to the peaks. These files are named by following the pattern of {Group_CD4CLAp0hct_vs_CD4CLAn0hct.DAR.txt}, in which CD4CLAp0hct anc CD4CLAn0hct are replaced by different groups existing in the covariate factor.


