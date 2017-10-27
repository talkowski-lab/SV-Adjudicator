# 04. Variant resolution
This module merges structural variation calls across all PE/SR and depth algorithms after filtering, then resolves complex variants from the consensus inversion/translocation breakpoints and CNV intervals.
## Process through snakemake

### Module configuration
The configuration file `config.yaml` outlines the module's inputs and parameters, and should be modified accordingly to each specific project. 

* `batches` : Sample/group/batch key.

* `samples` : list of individules to be included

* `chromos` : list of chromosomes to be processed.
This file should be modified according to different reference genome. It is recommended that autosomes and allosomes are prepared differently.

* `pesr_sources` : 	Names of pair end/split read algorithms to be processed

* `depth_sources` :	Names of read depth algorithms to be processed

* `input_vcfs` : folder containing vcf files from `03_variant_filtering`
../03_variant_filtering/filtered_vcfs/

* `pe_counts` : pe_counts/{batch}.pe.sorted.txt.gz 

* `sr_counts` : sr_counts/{batch}.sr.sorted.txt.gz

* `coveragefile`: {batch}.binCov.bed.gz 

* `medianfile` : {batch}.binCov.median

* `famfile` : ../../ref/{batch}.fam

* `cutoffs`: cutoff information trained from the Random Forest

### Input

### Output



