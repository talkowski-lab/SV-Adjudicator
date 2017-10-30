# PeTest
This repository contains the workflow that evaluates pair-end support for all SV calls on a per-batch basis.

## Required matrics
Pe matrics should be prepared for this process. The matircs describes discordant read pairs in all individuals, and can be collectd from the aligned sequenes by following these steps:

1. run `svtools collect-pesr` to collect discordant pair end information:	
```
svtools collect-pesr sample.bam split_count/sample.txt pe_count/sample.txt
```

2. add sample name as an extra column to each pe_count output:
```
python script/add_sample_name.py pe_count/sample.txt
```

Now the pe_count/sample.txt looks like this:
```
1	9995	-	1	249240383	-	sample
1	9997	-	3	197900294	-	sample
1	9999	-	1	10179	+	sample
1	9999	-	5	11254	+	sample
```
3. concatinate the pe_count output of all samples, sort, bgzip and tabix:
```
cat pe_count/*.txt > matircs.pe.txt
sort -k1,1 -k2,2n  matircs.pe.txt >  matircs.pe.sorted.txt
bgzip -c matircs.pe.sorted.txt > matircs.pe.sorted.txt.gz
tabix -b 2 -e 2 matircs.pe.sorted.txt.gz
```

## Input files
Input files for PeTest are produced through the 01_algorithm_integration step, and are usually kept under `../../01_algorithm_integration/vcfcluster/`. Names of input files are in this format: `{batch}.{source}.{chrom}.vcf.gz`

## Process through snakemake
You can simply modify the `config.yaml` file to fit your data, and type `snakemake` under this directory to get the RdTest workflow autonamously processed. Details as how to modify the `config.yaml` is described in **Module configuration and input**.

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

* `pe_counts` : pe_counts/{batch}.pe.sorted.txt.gz
This matrices contains all disconcordant pair end information are kept. To prepare this matrics, refer to **Required matrics** for instructions.

* `famfile` : ../../ref/{batch}.fam
This file describes the family structure in batch



## Process each script manually
Autosomes and allosomes should be processed separately, with two whitelists contaning samples names of all males (whitelists/{batch}.males.list) and females(whitelists/{batch}.females.list) prepared. The whitelists have one sample name in each line. 

For autosomes:
```
svtools pe-test ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz matircs.pe.sorted.txt.gz petest/{batch}.{source}.{chrom}.stats
```
For allosomes:
```
svtools pe-test --samples whitelists/{batch}.females.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.females.stats
svtools pe-test --samples whitelists/{batch}.males.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.males.stats
python script/pe_merge_allosomes.py batch source chrom X
python script/pe_merge_allosomes.py batch source chrom Y
```

## Efficient manual process
It is recommended that big vcf input be split randomly into smaller files (e.g. 100 SV records per vcf), run through pe-test and then merge the splits together. To split vcf files:
```
python script/split_vcf.random.py ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz split_vcf/{batch}.{source}.{chrom}. -s number_of_svs_per_split
```
Sort and index each split vcf (i.e. `split_out`),
```
vcf-sort split_out > split_out.vcf
bgzip split_out.vcf
tabix split_out.vcf.gz
```
For each split vcf, apply `svtools pe-test` and them merge them:
```
svtools pe-test split_out.vcf.gz matircs.pe.sorted.txt.gz split_petest/split_out.stats
cat {input} | sed -r -e '/^chr\\s/d' | sort -k1,1V -k2,2n | cat <(head -n1 {input[0]}) - > {output}
```

Here's full instruction of `script/split_vcf.random.py`:
```
usage: split_pesrtest_random.py [-h] [-s SIZE] input output

positional arguments:
  input                 namd of input vcf.gz to be splited
  output                prefix of output

optional arguments:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  size of outputs
```

## Output filesthi sformat: 
Result from this step are kept under the `petest/` folder with names in the format: `{batch}.{source}.{chrom}.stats` 



