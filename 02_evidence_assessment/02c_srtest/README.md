# SrTest
This repository contains the workflow that evaluates split read support for all SV calls on a per-batch basis.

## Required matrics
Sr matrics should be prepared for this process. The matircs describes soft clipped alignments in all individuals, and can be collectd from the aligned sequenes by following these steps:

1. run `svtools collect-pesr` to collect split read information:
```
svtools collect-pesr sample.bam split_count/sample.txt pe_count/sample.txt
```

2. add sample name as an extra column to each pe_count output:
```
python script/add_sample_name.py split_count/sample.txt
```

Now the split_count/sample.txt looks like this:
```
1	9997	left	1	sample
1	9999	left	3	sample
1	10000	left	1	sample
1	10001	left	5	sample
1	10002	left	9	sample
1	10003	left	10	sample
```

3. concatinate the split_count output of all samples, sort, bgzip and tabix:
```
cat split_count/*.txt > matircs.sr.txt
sort -k1,1 -k2,2n  matircs.sr.txt >  matircs.sr.sorted.txt
bgzip -c matircs.sr.sorted.txt > matircs.sr.sorted.txt.gz
tabix -b 2 -e 2 matircs.sr.sorted.txt.gz
```

## Input files
Input files for PeTest are produced through the 01_algorithm_integration step, and are usually kept under `../../01_algorithm_integration/vcfcluster/`. Names of input files are in this format: `{batch}.{source}.{chrom}.vcf.gz`

## Process through snakemake
You can simply modify the `config.yaml` file to fit your data, and type `snakemake` under this directory to get the RdTest workflow autonamously processed.

## Module configuration and input
The configuration file `config.yaml` outlines the module's inputs and parameters, and should be modified accordingly to each specific project. 

* `batches` : filepath
Sample/group/batch key.

* `input_vcfs` : vcf files to be processed in this step. 
The vcf files are produced through `01_algorithm_integration` and are kept under `../../01_algorithm_integration/vcfcluster`

* `input_beds` : vcf files to be processed in this step. 
The bed files are produced through `01_algorithm_integration` and are kept under `../../01_algorithm_integration/rdtest_beds`

* `groups` : list of samples to be processed.

* `chromos` : list of chromosomes to be processed.
This file should be modified according to different reference genome. It is recommended that autosomes and allosomes are prepared differently.

* `pesr_sources` : 
Names of pair end/split read algorithms to be processed

* `depth_sources` :
Names of read depth algorithms to be processed

*`sr_counts` : pe_counts/{batch}.sr.sorted.txt.gz
This matrices contains all split read information. To prepare this matrics, refer to **Required matrics** for instructions.

* `famfile` : ../../ref/{batch}.fam
This file describes the family structure in batch

## Process each script manually
Autosomes and allosomes should be processed separately, with two whitelists contaning samples names of all males (whitelists/{batch}.males.list) and females(whitelists/{batch}.females.list) prepared. The whitelists have one sample name in each line.
For autosomes:
```
svtools pe-test ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz matircs.sr.sorted.txt.gz srtest/{batch}.{source}.{chrom}.stats
```
For allosomes:
```
svtools sr-test --samples whitelists/{batch}.females.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.sr.sorted.txt.gz srtest_allosomes/{batch}.{source}.{chrom}.females.stats
svtools sr-test --samples whitelists/{batch}.males.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.sr.sorted.txt.gz srtest_allosomes/{batch}.{source}.{chrom}.males.stats
python script/sr_merge_allosomes.py batch source chrom X
python script/sr_merge_allosomes.py batch source chrom Y
```

## Output filesthi sformat:
Result from this step are kept under the `srtest/` folder with names in the format: `{batch}.{source}.{chrom}.stats`

