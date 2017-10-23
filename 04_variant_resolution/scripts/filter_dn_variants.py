#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import pandas as pd
import pysam
import svtools.utils as svu


def set_null(record, sample):
    dat = record.samples[sample].items()

    # Set genotype to no-call
    record.samples[sample]['GT'] = (0, 0)

    for fmt, value in dat:
        if fmt == 'GT':
            continue

        # Get type and count of FORMAT
        n = record.format[fmt].number
        dtype = record.format[fmt].type

        # Set null value based on FORMAT type
        if dtype in 'Integer Float'.split():
            null_val = 0
        elif dtype in 'String Character'.split():
            null_val = ''
        else:
            raise ValueError('Invalid VCF FORMAT type: {0}'.format(dtype))

        # Set the appropriate count of values
        if n == 1:
            record.samples[sample][fmt] = null_val
        elif n == '.':
            record.samples[sample][fmt] = (null_val, )
        else:
            record.samples[sample][fmt] = tuple(null_val for i in range(n))


def filter_dn_variants(vcf, metrics, fout):
    CNV = 'DEL DUP'.split()
    is_cnv = metrics.svtype.isin(CNV)
    depth_only = metrics.sources == 'depth'
    pesr_size_filter = (metrics.svsize >= 1000)

    passing = ((depth_only & is_cnv & metrics.rd_pass) |
               (~depth_only & is_cnv & pesr_size_filter & metrics.rd_pass) |
               (~depth_only & is_cnv & ~pesr_size_filter & metrics.pesr_pass) |
               (~is_cnv & metrics.pesr_pass))

    def _join_samples(s):
        return sorted(set(s))

    passes = metrics.loc[passing].groupby('name')['sample'].agg(_join_samples)
    fails = metrics.loc[~passing].groupby('name')['sample'].agg(_join_samples)

    checked_variants = metrics.name.unique()
    for record in vcf:
        # Write records unaltered if they weren't included in the de novo check
        if record.id not in checked_variants:
            fout.write(record)
        # Otherwise set samples appropriately
        else:
            pass_samples = passes.get(record.id, [])
            for sample in pass_samples:
                record.samples[sample]['GT'] = (0, 1)

            fail_samples = fails.get(record.id, [])
            for sample in fail_samples:
                set_null(record, sample)

            # Only report record if any samples made it through de novo check
            if len(svu.get_called_samples(record)) > 0:
                fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('metrics', help='Classified metrics')
    parser.add_argument('fout', help='Filtered VCF')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    fout = sys.stdout if args.fout in 'stdout -'.split() else args.fout
    fout = pysam.VariantFile(fout, 'w', header=vcf.header)

    metrics = pd.read_table(args.metrics)

    filter_dn_variants(vcf, metrics, fout)


if __name__ == '__main__':
    main()
