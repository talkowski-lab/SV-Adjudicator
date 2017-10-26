## Process through snakemake
Just type `snakemake` under this directory and the RdTest workflow with be autonamously processed.

## Process each step manually
It is also possible to run each step in this module manually, by following these steps:

1. Randomly split the big bed file into small beds for faster processing:
```
python scripts/split_rdtest_random.py ../../01_algorithm_integration/rdtest_beds/{batch}.{source}.{chrom}.bed split_beds/{batch}.{source}.{chrom}. -s number_SVs_per_bed
```
Output from this step are `{batch}.{source}.{chrom}.XXX`, where XXX are numbers differentiating each sub-file.

2. Run RdTest. In this step, samples should be split by sex to have chromosome X and Y processed. 
For **autosomes**:

```
Rscript scripts/RdTest.R -b split_beds/{batch}.{source}.{chrom}.XXX -o split_rdtest/ -n {batch}.{source}.{chrom}.XXX -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam
```

For **allosomes**, a female whitelist containing full list of all females samples and a male list with names of all males samples should be prepared:
```
Rscript ./scripts/RdTest.R -b split_beds/{batch}.{source}.{chrom}.XXX -o split_rdtest/ -n {batch}.{source}.{chrom}.XXX.females -w female.whitelist  -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam
Rscript ./scripts/RdTest.R -b split_beds/{batch}.{source}.{chrom}.XXX -o split_rdtest/ -n {batch}.{source}.{chrom}.XXX.males -w male.whitelist  -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam
```

3. Concatinate the sub metrics:
```
bash scripts/rdtest_mergesplit.sh batch source chrom
```


## Output files
Output files include RdTest scores for each CNV, are kept under `rdtest/`:
* `rdtest/{batch}.{source}.{chrom}.metrics`

