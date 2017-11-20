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

def remove_outliers(vcf, outliers):
    """
    Remove records specific to outlier samples

    Parameters
    ----------
    vcf : pysam.VariantFile
        VCF to filter
    outliers : dict of {svtype (str): samples (list of str)}
        List of outlier sample IDs for each svtype

    Yields
    ------
    record : pysam.VariantRecord
    """

    # TODO: skip filtering if no samples outliers
    # TODO: skip filtering per svtype if no samples outliers
    samples = list(vcf.header.samples)
    null_GTs = [(0, 0), (0, ), (None, None), (None, )]

    for record in vcf:
        svtype = record.info['SVTYPE']

        # Yield if no outliers specified
        if svtype not in outliers.keys():
            yield record
            continue

        # If no outliers present, return record
        outlier_samples = [s for s in samples if s in outliers[svtype]]
        if len(outlier_samples) == 0:
            yield record
            continue

        # Check if any normal samples are called
        def _is_called(sample):
            return record.samples[sample]['GT'] not in null_GTs
        normal_samples = [s for s in samples if s not in outliers[svtype]]
        normal_called = any([_is_called(s) for s in normal_samples])

        # Exclude records called only in outliers or not called in any samples
        if normal_called:
            yield record


def main():
    parser = argparse.ArgumentParser("find_outliers.py")
    parser.add_argument("input", type=str, help="list of samples names")
    parser.add_argument("output", type=str, help="list of samples names")
    parser.add_argument("outliers", type=str, help="list of samples names")

    args = parser.parse_args()

    #vcf = VariantFile(snakemake.input.vcf)
    vcf = VariantFile(args.input)

    outlier_table = pd.read_table(args.outliers)
    filtered = VariantFile(args.output, mode='w', header=vcf.header)

    outliers = defaultdict(list)
    for idx, row in outlier_table.iterrows():
        outliers[row['svtype']].append(row['sample'])

    for record in remove_outliers(vcf, outliers):
        filtered.write(record)

    filtered.close()


if __name__ == '__main__':
    main()

