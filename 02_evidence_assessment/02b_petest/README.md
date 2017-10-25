#PeTest
This repository contains the workflow that evaluates pair-end support for all SV calls on a per-batch basis.

##Required matrics
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

##Input files
Input files for PeTest are produced through the 01_algorithm_integration step, and are usually kept under `../../01_algorithm_integration/vcfcluster/`. Names of input files are in this format: `{batch}.{source}.{chrom}.vcf.gz`


##Modify config.yaml

##Run through snakemake
Just type snakemake under this directory and the RdTest workflow with be autonamously processed.

##Run each script manually
Autosomes and allosomes should be processed separately, with two whitelists contaning samples names of all males (whitelists/{batch}.males.list) and females(whitelists/{batch}.females.list) prepared. The whitelists have one sample name in each line. 

For autosomes:
```
svtools pe-test ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz matircs.pe.sorted.txt.gz petest/{batch}.{source}.{chrom}.stats

```
For allosomes:
```
svtools pe-test --samples whitelists/{batch}.females.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.females.stats
svtools pe-test --samples whitelists/{batch}.males.list ../../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz  matircs.pe.sorted.txt.gz petest_allosomes/{batch}.{source}.{chrom}.males.stats
python script/merge_allosomes.py batch source chrom X
python script/merge_allosomes.py batch source chrom Y
```

##Output filesthi sformat: 
Result from this step are kept under the `petest/` folder with names in the format: `{batch}.{source}.{chrom}.stats` 



