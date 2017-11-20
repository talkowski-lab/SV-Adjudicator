#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
find_outliers.py

Identify per-algorithm, per-svtype outliers
"""


from pysam import VariantFile
from collections import defaultdict
import pandas as pd
import os
import argparse
def find_outliers(vcflist, svtypes):
    """
    Locate samples which have abnormally high count of a specific variant type.

    Parameters
    ----------
    vcflist : list of str
        Filepaths to VCFs
    svtype : str
        SV class to count [DEL,DUP,INV,BND]

    Returns
    -------
    outliers : list of str
        List of outlier samples
    """

    counts = defaultdict(int)
    null_GTs = [(0, 0), (0, ), (None, None), (None, )]

    # Count calls per sample of specified svtype
    for f in vcflist:
        vcf = VariantFile(f)
        for record in vcf:
            svtype = record.info['SVTYPE']
            if svtype not in svtypes:
                continue

            for sample in record.samples:
                gt = record.samples[sample]['GT']
                if gt not in null_GTs:
                    counts[(sample, svtype)] += 1

    counts = pd.DataFrame.from_dict({'var_count': counts})
    columns = {'level_0': 'sample', 'level_1': 'svtype'}
    counts = counts.reset_index().rename(columns=columns)

    # Calculate outlier cutoff
    Q1 = counts.groupby('svtype').var_count.quantile(0.25)
    Q3 = counts.groupby('svtype').var_count.quantile(0.75)
    IQR = Q3 - Q1
    cutoff = Q3 + 1.5 * IQR

    # Identify outliers
    cutoff = cutoff.rename('cutoff').reset_index()
    outliers = pd.merge(counts, cutoff, on='svtype', how='left')
    outliers = outliers.loc[outliers.var_count > outliers.cutoff]

    return outliers


def main():
    parser = argparse.ArgumentParser("find_outliers.py")
    parser.add_argument("input", type=str, help="list of samples names")
    parser.add_argument("output", type=str, help="list of samples names")
    parser.add_argument("svtype", type=str, help="list of samples names")

    args = parser.parse_args()

    vcfprefix = args.input.split('/')[-1]
    vcflist=[]
    for i in os.listdir('std_vcfs/'):
          if vcfprefix in i:   vcflist.append('std_vcfs/'+i)

    outliers = find_outliers(vcflist, args.svtype.split(','))
    outliers.to_csv(args.output, index=False, sep='\t')


if __name__ == '__main__':
    main()

