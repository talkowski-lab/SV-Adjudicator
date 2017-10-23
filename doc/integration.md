# Algorithm integration

The algorithm integration module combines variant predictions across each set
of algorithms (i.e. PE/SR and RD), then integrates across the two classes of
evidence.

## PE/SR integration
All data based on the standardized calls provided in 
`{workdir}/preprocessing/filtered_vcfs/{source}.{quad}.vcf.gz`

* `vcflists/vcfs.{source}.list`  
    List of VCFs being clustered
* `vcfcluster/{source}/merged.{chrom}.vcf`  
    Clustered variants.

## RD integration
All data based on the standardized calls provided in 
`{workdir}/preprocessing/std_beds/merged.{svtype}.bed.gz`

Per-file column details to be described later as necessary.

* `depth_intersect/merged.{svtype}.{chrom}.bed.gz`  
    Self-intersection of standardized calls, split by chromsoome. 
* `depth_variant_lists/merged.{svtype}.{chrom}.list`  
    List of variant IDs present in standardized calls, split by chromosome.
    Necessary for bedcluster's sparse graph clustering.
* `bedcluster/merged.{chrom}.svof`  
    Clustering of variants, linking variants based on self-intersection hits

## RdTest filtering prior to PE/SR+RD integration

* `rdtest_beds/{source}/merged.{chrom}.bed`  
    Variants formatted for RdTest. 
* `rdtest_filtered/depth/merged.{chrom}.svof`  
    Clustered depth calls which passed RdTest thresholds.
