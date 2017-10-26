# 03 Variant Filtering
This workflow integrate, filters and genotypes the structural variation(SVs) calls based on the evidence derived in previous modules. The following processes are applied here:
1. Evidences collected from module 02, i.e., rd, pe, sr and baf, are aggregated;
2. Process integrated evidences through Random Forest(RF) to train the optimized parameter for quality control
3. Apply the RF filters on the SVs and remove those failures.

## Process through snakemake
### Module configuration
### Input
VCFs and BEDs integrated from `01_algorithm_integration` and evidences (i.e. rd, pe, sr and baf) matrices produced from `02_evidence_assessment` are required for this workflow:
* rdtest: ../02_evidence_assessment/02a_rdtest/rdtest
* petest: ../02_evidence_assessment/02b_petest/petest
* srtest: ../02_evidence_assessment/02c_srtest/srtest
* baftest: ../02_evidence_assessment/02d_baftest/baftest
* input_vcfs: ../01_algorithm_integration/vcfcluster
* input_beds: ../01_algorithm_integration/rdtest_beds

A batch_list including names of all samples in the same batch should be prepared for read depth calls to be properly aggregated
* batch_list

## Manual process
#### Evidence aggragation
a. To aggregate evidence for **pesr callers** (eg. delly, lumpy, manta, wham), for each `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-r ../02_evidence_assessment/02a_rdtest/rdtest/{batch}.{source}.{chrom}.metrics \
	-p ../02_evidence_assessment/02b_petest/petest/{batch}.{source}.{chrom}.stats \
	-s ../02_evidence_assessment/02c_srtest/srtest/{batch}.{source}.{chrom}.stats \
	-b ../02_evidence_assessment/02d_baftest/baftest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz \
	metrics/{batch}.{source}.{chrom}.metrics
```

b. To aggregate evidence for **read depth callers** (eg. cn.mops, CNVnator, ERDs), for `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-r ../02_evidence_assessment/02a_rdtest/rdtest/{batch}.{source}.{chrom}.metrics \
	-b ../02_evidence_assessment/02d_baftest/baftest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/rdtest_beds/{batch}.{source}.{chrom}.bed \
	--bed \
	metrics/{batch}.{source}.{chrom}.metrics
```

c. To aggregate evidence for **mobile element insertion callers** (eg. MELT), for `{source}` and `{chrom}`: 
```
python scripts/aggregate.py \
	-s ../02_evidence_assessment/02c_srtest/srtest/{batch}.{source}.{chrom}.stats \
	-v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz \
	--batch-list batch_list
	metrics/{batch}.{source}.{chrom}.metrics
```

Modify the aggragated metrics with position of each variants
```
python scripts/add_pos.py -v ../01_algorithm_integration/vcfcluster/{batch}.{source}.{chrom}.vcf.gz metrics/{batch}.{source}.{chrom}.metrics
python scripts/add_pos.py -v ../01_algorithm_integration/rdtest_beds/{batch}.{source}.{chrom}.bed metrics/{batch}.{source}.{chrom}.metrics
