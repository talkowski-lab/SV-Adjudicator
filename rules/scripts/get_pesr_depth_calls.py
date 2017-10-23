#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import pandas as pd


def group_samples(df, source):
    """
    Make a list of samples called in each variant.
    """
    samples = df.groupby('name')['sample']\
                .agg(lambda s: sorted(set(s)))\
                .rename('{0}_samples'.format(source)).reset_index()\
                .rename(columns={'name': '{0}_name'.format(source)})

    return samples


def union_samples(row):
    """Take union of samples"""
    pesr = row.pesr_samples
    depth = row.depth_samples
    samples = sorted(set(pesr + depth))
    return ','.join(samples)


def get_pesr_depth_calls(pesr, depth, links):
    """
    Take the union of samples in linked PE/SR and depth cals
    """

    # Make unique list of variants based on PE/SR coordinates and
    # add linked depth names
    variants = pesr['#chrom start end name svtype'.split()].drop_duplicates()
    variants = variants.loc[variants.name.isin(links.pesr_name)].copy()
    variants['pesr_name'] = variants['name']
    variants = pd.merge(variants, links, on='pesr_name', how='left')

    # Add PE/SR and depth samples to each variant
    pesr_samples = group_samples(pesr, 'pesr')
    depth_samples = group_samples(depth, 'depth')
    variants = pd.merge(variants, pesr_samples, on='pesr_name', how='left')
    variants = pd.merge(variants, depth_samples, on='depth_name', how='left')

    variants['samples'] = variants.apply(union_samples, axis=1)

    cols = '#chrom start end name svtype samples pesr_name depth_name'.split()
    return variants[cols]


def main():
    cols = '#chrom start end name svtype sample'.split()
    pesr = pd.read_table(snakemake.input.pesr)
    depth = pd.read_table(snakemake.input.depth)[cols].copy()
    links = pd.read_table(snakemake.input.links,
                          names='pesr_name depth_name'.split())

    pesr_depth = get_pesr_depth_calls(pesr, depth, links)
    pesr_depth.to_csv(snakemake.output[0], sep='\t', index=False)


if __name__ == '__main__':
    main()
