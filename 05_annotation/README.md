# 05. Variant annotation
This module annotate the integrated and filtered variants ... ...

## Process through snakemake

### Module configuration
The configuration file `config.yaml` outlines the module's inputs and parameters, and should be modified accordingly to each specific project. 

* `gencode` : 

* `annotation_gtf` :

* `pc_transcripts_fa` :

* `pc_translations_fa` :

* `transcript_source` :

Annotation files for hg19 are provided under this folder `annotation` as an example. For other versions of reference, correspondance annotations are availble at ... ...

### Input
Result from previous step should be provided as input for this step.

### Output
Annotated variants are reported in vcf format from this workflow.

## Manual process
1a. produce gencode.canonical_annotation

```
python scripts/get_canonical_transcripts.py \
		annotation/annotation_gtf \
		annotation/pc_translations_fa \
		annotation/pc_transcripts_fa \
		annotation/transcript_source \
		gencode/gencode.canonical_transcripts.txt

cat 
    <(cut -f2 gencode/gencode.canonical_transcripts.txt | sed -e '1d' | fgrep -w -f - <(zcat annotation/annotation_gtf)) \
    <(cut -f1 gencode/gencode.canonical_transcripts.txt | sed -e '1d' | fgrep -w -f - <(zcat annotation/annotation_gtf) | awk '($3=="gene")') \
  | sort -k1,1V -k4,4n \
  | sed -e 's/^chr//' \
  | bgzip -c \
  > gencode/gencode.canonical_annotation.gtf.gz

```

1b. produce 
... ...

2. annotate:
```
svtype annotate \
	--gencode gencode/gencode.canonical_annotation.gtf.gz \
	--noncoding noncoding/noncoding_elements.bed \
	vcfs/variants.{chrom}.vcf.gz \
	annotated_vcfs/variants.{chrom}.vcf
```

