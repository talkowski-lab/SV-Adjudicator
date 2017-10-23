#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import os
from collections import deque
import pandas as pd


def get_coord(row):
    if row.svtype == 'DEL':
        if row['clip'] == 'right':
            return 'start'
        elif row['clip'] == 'left':
            return 'end'
    elif row.svtype == 'DUP':
        if row['clip'] == 'right':
            return 'end'
        elif row['clip'] == 'left':
            return 'start'


def coord_dist(row):
    coord = row['coord']

    variant_pos = row[coord]
    return variant_pos - row.pos


def fill_missing_dists(splits, samples):
    dists = splits['dist'].unique()
    idx = pd.MultiIndex.from_product(iterables=[samples, dists])

    splits = splits.set_index('sample dist'.split())
    splits = splits.reindex(idx).fillna(0).astype(int).reset_index()
    splits = splits.rename(columns=dict(level_0='sample', level_1='dist'))

    return splits


def main():
    # get list of samples called in each variant

    # Read in each dataframe (until pysam can read from S3)
    bed = pd.read_table(snakemake.input.bed).drop('samples', axis=1)

    sample_list = pd.read_table(snakemake.input.samples)
    samples = sample_list['sample'].unique()

    names = 'pos count clip sample'.split()
    splits = pd.read_table(snakemake.input.counts, names=names)
    splits = splits.drop_duplicates()
    splits = splits.loc[splits['sample'].isin(samples)].copy()

    splits['name'] = snakemake.wildcards.name
    splits = pd.merge(splits, bed, on='name', how='left')

    # Calculate distance from start/end
    splits['svtype'] = splits.svtype.str.upper()
    splits['coord'] = splits.apply(get_coord, axis=1)
    splits['dist'] = splits.apply(coord_dist, axis=1)

    # Filter by distance to appropriate coordinate.
    # Accounts for cases where a split is within the window around one
    # coordinate, but is clipped in the opposite direction
    window = snakemake.params.window
    splits = splits.loc[splits['dist'].abs() < window].copy()

    cols = 'sample coord dist count'.split()
    splits = splits[cols].copy()

    # Fill in missing samples at any pos with a split read in a sample
    sub_splits = []
    for coord in 'start end'.split():
        df = splits.loc[splits.coord == coord, 'sample dist count'.split()]
        df = fill_missing_dists(df, samples)
        df['coord'] = coord
        sub_splits.append(df)
    splits = pd.concat(sub_splits)

    cols = 'sample call_status'.split()
    splits = pd.merge(splits, sample_list[cols], on='sample', how='left')

    splits['name'] = snakemake.wildcards.name
    splits.to_csv(snakemake.output[0], index=False, sep='\t')


if __name__ == '__main__':
    main()
