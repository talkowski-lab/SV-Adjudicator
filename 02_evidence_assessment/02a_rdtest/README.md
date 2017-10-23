# RdTest

This workflow evaluates read-depth support for all CNV calls on a per-batch
basis.

## Input files

* `input_beds/{batch}.{source}.{chrom}.bed`
    CNV calls in RdTest format.  
    `#chrom	start	end	name	samples	svtype`

## Output files (in progress)

* `rdtest/{batch}.{source}.{chrom}.bed.pk`
    RdTest scores for each CNV. For further details see RdTest documentation.
