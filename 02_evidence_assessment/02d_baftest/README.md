# BAFTest
This repository contains the workflow to generates **BAF reference file **, i.e. gzipped & tabix indexed file containing all BAF information needed for BAF analysis, and process the BAF test of all variants

## Required matrics
BAF metrics should be prepared. Follow theses steps for the generation: 

### Required Input:
   1. A gzipped & tabix indexed vcf file (s3 bucket compatible)
   2. A fasta.fai index file (can be the hg19 fai file, or if needed, select chromosomes), to split the genome into manageble chunks

### Step 1. Chunk the genome
 ```
 python ./script/splitchr.py ${refdict} 500000 > chunk.txt
 ```

### Step 2. Extract BAF information from the chunks
```
touch mytest1.snp
tabix -h ${vcf} ${chr}:${start}-${end} | grep -E "^#|PASS"  |bcftools view -M2 -v snps - |python scratch.py
 ```

### Step 3. Combine the BAF info chunks, gzip and index
 ```cat ${sep=" " bafs} > baf_snp.txt
    sort -k1,1d -k2,2n baf_snp.txt>baf_snp_sorted.txt
    bgzip baf_snp_sorted.txt
    tabix -b2 baf_snp_sorted.txt.gz
 ```


#### NOTE:
 Batch information & ID-swapping are currently disabled. Assuming the IDs stay consistent throughout analysis


## Apply BAF analysis
### Input: 
baf file generated in the previous step (baf_snp_sorted.txt.gz), a standarized bed file listing the variants to validate

### Example
 ```
 python scrat.py Phase1.delly.22.stat ../Filegenerate/baf_snp_sorted.txt.gz --batch Phase1 > baf.result
 ```


