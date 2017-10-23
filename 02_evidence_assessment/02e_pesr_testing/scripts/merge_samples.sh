#! /bin/bash
#
# merge_samples.sh
# Copyright (C) 2017 ec2-user <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.
#


samples=$1
chrom=$2

cmd="sort -m -k1,1V -k2,2n -k4,4V -k5,5n"

while read sample; do
  cmd="$cmd <(bgzip -d -c sample_counts/${chrom}/${sample}.txt.gz)"
done < $samples

echo $cmd

mkdir -p merged_counts
eval "$cmd" | bgzip -c > merged_counts/cohort.${chrom}.txt.gz
