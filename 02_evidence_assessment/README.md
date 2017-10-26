# 02. Evidence assessment

This module provides workflows for assessing read-depth (RD), paired-end (PE),
split read (SR), and B-allele frequency (BAF) evidence supporting strutural
variant calls.

## Required matrics
To assess the evidence for potential SVs, there are three matrics requred as input:
* `Matrics.binCov.bed.gz` and `Matrics.binCov.median`
* `Matrics.pe.sorted.txt.gz`
* `Matrics.sr.sorted.txt.gz`

All matrics should have been properlly bgziped and tabix indexed. Refer to each subdirectory for details of how these matrics are created and applied to the analysis.

## Quick process through snakemake
Each evidence can be processed individually by running *snakemake* under each sub-directory, e.g.:
```
cd 02a_rdtest
snakemake
```

For futher details for processing each evidence, please refer to the readme under each sub-directory.
