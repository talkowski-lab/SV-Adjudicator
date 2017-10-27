# Preprocessing

This workflow standardizes SV calls across a series of algorithms to a
consistent VCF format for PE/SR algorithms and a consistent BED format for
depth algorithms.

This workflow contains many elements which are specific to the SSC analyses for
which the pipeline was initially developed, and is not intended to be
generalizable to all use cases. However, it is presented here for a)
reproducibility of the SSC analyses and b) as an example of how SV calls may be
standardized to the format required by the remaining steps of the pipeline.

## Process through snakemake
### Module configuration and input
The configuration file `config.yaml` outlines the module's inputs and parameters, and should be modified accordingly to each specific project. 

* `batches` : filepath
Sample/group/batch key.

* `groups` : filepath
List of groups to include during integration. Expects one PE/SR VCF per group.

* `pesr_sources` : 
Names of pair end/split read algorithms to be processed

* `depth_sources` :
Names of read depth algorithms to be processed

* `svtypes` :
Types of SVs to be processed

* `cnv_types` :
Types of CNVs to be processed

* `outlier_removal` :
Sources to remove outliers from 

### Input data

#### PE/SR calls
PE/SR calls should be placed in the `data/raw_vcfs/` directory and follow the
filename convention `{source}/{source}.{group}.vcf.gz`, where `source` refers to the
source algorithm which generated the calls and `group` refers to an identifier
of the group of samples on which the algorithm was run. 

#### RD calls
RD calls across all samples should be concatinated together, and placed in tbe `data/raw_beds/` directory and follow the
filename convention `{source}/{sample}.{svtype}.raw.bed`. `svtype` can be either DEL or DUP.

### Output data

#### PE/SR preprocessing
The PE/SR preprocessing module standardizes VCFs and removes calls specific to
outlying samples. Outliers are determined to be samples with more than 
`(Q3 + 1.5 * IQR)` variants observed, and are calculated on a per-algorithm,
per-svtype basis (e.g. Delly inversions).

* `std_vcfs/{source}.{quad}.vcf`  
    Standardized VCFs.  
        - INFO fields for CHR2, END, STRANDS, SVTYPE, SVLEN, and SOURCE    
        - BND ALTs are converted to VCF specification  
        - Balanced events are separated into stranded breakpoints  

* `outliers/{source}.list`  
    Tables of svtype-specific outlier samples by caller.  
        - `sample`: sample ID  
        - `svtype`: type of SV for which the sample is an outlier  
        - `var_count`: number of variants observed in the sample  
        - `cutoff`: Q3 + 1.5 * IQR for that svtype/algorithm  

* `filtered_vcfs/{source}.{quad}.vcf.gz`  
    VCFs with calls specific to outlier samples or null in all samples removed.
    Sorted and tabix-indexed to prepare for clustering.

#### RD preprocessing
The RD preprocessing module merges all calls made by any depth algorithm in a 
sample to create a single bed per sample, annotated with the algorithms which
contributed to each merged call.

* `std_beds/`
    * `{source}/{sample}.{svtype}.raw.bed`  
        Raw calls. Columns are standardized to chrom, start, end, svsize
    * `{source}/{sample}.{svtype}.merged.bed`  
        Raw beds subjected to bedtools merged.
    * `merged_algs/{sample}.{svtype}.raw.bed`  
        Sorted concatenation of raw calls.
    * `merged_algs/{sample}.{svtype}.merged.bed`  
        Per-sample merging of raw calls across all algorithms.
    * `merged_algs/{sample}.{svtype}.cov.bed`  
        Final per-sample standardized bed. Each merged call is queried for
        coverage by each individual algorithm's callset in order to determine 
        the algorithms which predicted the variant.  Columns are chrom, start, 
        end, sample svsize, cnmops coverage, cnvnator coverage, and genomestrip 
        coverage. 
    * `{batch}.{svtype}.bed.gz`  
        Concatenation of standardized `*.cov.bed` from all samples in a batch.

## Manual process
### PE/SR
PE/SR preprocess includes these steps: 
1. standardize vcf format
```
svtools standardize  {source}.{group} data/raw_vcfs/{source}/{source}.{group}.vcf.gz std_vcfs/{source}.{group}.vcf {source}
```
2. identify and remove vcf outliers  
```
python scripts/find_outliers.py std_vcfs/{source}.{group}.vcf outliers/{source}.list
python scripts/remove_outliers.py std_vcfs/{source}.{group}.vcf outliers/{source}.list filtered_vcfs/{source}.{group}.vcf
```
3. sort and index vcf
```
vcf-sort -c filtered_vcfs/{source}.{group}.vcf | bgzip -c > filtered_vcfs/{source}.{group}.vcf.gz
tabix -p filtered_vcfs/{source}.{group}.vcf.gz
```
### RD
RD preprocess includes:
1. standardize vcf files from individual algorithm:
cnmops:
```
zcat data/raw_beds/cnmops/{sample}.{svtype}.raw.bed \
	| fgrep {sample} \
	| fgrep {svtype} \
	| awk -v OFS="\\t" '($1!="X" && $1!="Y") {{print $1, $2, $3, $3-$2, "1", "0", "0"}}' \
	| cat - <(sed -e '/^#/d' {input.allosomal} | awk -v OFS="\\t" '{{print $1, $2, $3, $3-$2, 1, 0, 0}}') \
	> std_beds/cnmops/{sample}.{svtype}.raw.bed
```
cnvnator
```
awk -v OFS="\\t" '{{print $1, $2+50, $5+50, $8, "0", "1", "0"}}'  data/raw_beds/cnvnator/{sample}.{svtype}.win100.bp200.sort.bedpe > std_beds/cnvnator/{sample}.{svtype}.raw.bed
```
genomestrip
```
awk -v OFS="\\t" '{{print $1, $2, $3, $4-$3, "0", "0", "1"}}' data/raw_beds/genomestrip/{sample}.genomestrip.{svtype}.bed > std_beds/genomestrip/{sample}.{svtype}.raw.bed
```
2. concatinate standardized beds from each algorithm:
```
cat std_beds/cnmops/{sample}.{svtype}.raw.bed std_beds/cnvnator/{sample}.{svtype}.raw.bed std_beds/genomestrip/{sample}.{svtype}.raw.bed | sort -k1,1V -k2,2n > std_beds/merged_algs/{sample}.{svtype}.raw.bed
```

