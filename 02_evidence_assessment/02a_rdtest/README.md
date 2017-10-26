# RdTest

This repository contains the workflow that evaluates read-depth support for all CNV calls on a per-batch basis.

## Required matrics

Two matrics are required to process the RdTest: the Matrics.binCov.bed.gz and Matrics.binCov.median

* `Matrics.binCov.bed.gz` is a bgziped bed file that contains the bincov coverages accross whole genome of all individuals involved in the SV discovery project. This file should be tabix indexed. The first three columns describes the gennomic location and following columns describes the bincov coverage of each individual, e.g. 

|#chr | start | end | sample1 | sample2 | sample3 | sample4 | sample5 | sample6 | sample7 | sample8 | sample9 | sample10 | 
|-----|-------|-----|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|----------------------|
|1 | 10000 | 10100 | 938 | 1387 | 954 | 1344 | 688 | 1219 | 1662 | 2027 | 1221 | 1233|
|1 | 10100 | 10200 | 1089 | 1462 | 927 | 1365 | 840 | 1290 | 1774 | 2217 | 1316 | 1379|
|1 | 10200 | 10300 | 554 | 694 | 462 | 679 | 516 | 692 | 914 | 1140 | 787 | 637|
|1 | 10300 | 10400 | 1149 | 1473 | 1019 | 1458 | 945 | 1391 | 1977 | 2231 | 1433 | 1330|
|1 | 10400 | 10500 | 767 | 964 | 649 | 880 | 583 | 891 | 1336 | 1507 | 797 | 849|
|1 | 10500 | 10600 | 43 | 102 | 81 | 38 | 59 | 49 | 157 | 145 | 56 | 41|
|1 | 10600 | 10700 | 16 | 40 | 32 | 21 | 23 | 18 | 62 | 43 | 24 | 15|
|1 | 10700 | 10800 | 0 | 1 | 1 | 1 | 0 | 0 | 4 | 3 | 1 | 0|
|1 | 10800 | 10900 | 4 | 0 | 3 | 21 | 0 | 27 | 15 | 20 | 2 | 21|
|1 | 10900 | 11000 | 14 | 0 | 5 | 41 | 5 | 63 | 36 | 45 | 6 | 40|
|1 | 11000 | 11100 | 24 | 0 | 10 | 69 | 8 | 105 | 51 | 74 | 9 | 62|


* `Matrics.binCov.median`  contains the median of bincov coverages for each individual, e.g.

|#chr | start | end | sample1 | sample2 | sample3 | sample4 | sample5 | sample6 | sample7 | sample8 | sample9 | sample10 | 
|-----|-------|-----|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|----------------------|
|61 | 67 | 62 | 63 | 61 | 66 | 63 | 61 | 75 | 69 | 76 | 64 | 68 |


## Input files

input files for RdTest are produced through the *01_algorithm_integration* step, and are usually kept under: 
```
../../01_algorithm_integration/rdtest_beds/
```

names of input files are in this format: `{batch}.{source}.{chrom}.bed`, with the following example representing the actual information inside:
```
#chrom	start	end	name	samples	svtype
2	1	10268	CMC_depth_DEL_2_0	sample2   DEL
2	10000	18000	CMC_depth_DUP_2_1	sample8,sample10   DUP
2	10000	41072	CMC_depth_DUP_2_2	sample1,sample2,sample4,sample6,sample7,sample8,sample10
```

## Modify config.yaml


## Run through snakemake
Just type `snakemake` under this directory and the RdTest workflow with be autonamously processed.

## Run each script manually
It is also possible to run each step in this module manually, which allows user to process multiple jobs in parallel. 

First step is to randomly split the big bed file into small beds for faster processing:
```
python ./script/split_rdtest_random.py input.bed output.prefix -s number_SVs_per_bed
```

Output from this step are output.prefix.XXX, where XXX are numbers differentiating each sub-file. Next step is to run RdTest. In this step, autosomes and allosomes are treated differently, and output.prefix.XXX.metrics will be produced. 

For autosomes:

```
Rscript ./scripts/RdTest.R -b output.prefix.XXX -o output/path/ -n output.prefix.XXX  -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam
```

For allosomes, a female whitelist containing full list of all females samples and a male list with names of all males samples should be prepared:
```
Rscript ./scripts/RdTest.R -b output.prefix.XXX -o output/path/ -n output.prefix.XXX.females -w female.whitelist  -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam
Rscript ./scripts/RdTest.R -b output.prefix.XXX -o output/path/ -n output.prefix.XXX.males -w male.whitelist  -c Matrics.binCov.bed.gz -m Matrics.binCov.median -f full_sample_list.fam

another script to merge male and female calls to be added here .....
```

Last step is to concatinate the sub metrics:
```

```




## Output files (in progress)

* `rdtest/{batch}.{source}.{chrom}.bed.pk`
    RdTest scores for each CNV. For further details see RdTest documentation.
