#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import pandas as pd


def main():
    dfs = []
    for metricfile in snakemake.input.metrics:
        dfs.append(pd.read_table(metricfile))
    metrics = pd.concat(dfs)

    # Exclude variants where all samples are called
    model_mask = metrics.RD_Model != 'All_samples_called_CNV_no_analysis'
    metrics = metrics.loc[model_mask].copy()

    # Add CNV size and seg dup coverage
    coverage = pd.read_table(snakemake.input.coverage)
    metrics = pd.merge(metrics, coverage, on='name', how='left')

    # Sort by name
    metrics = metrics.sort_values('name')

    metrics.to_csv(snakemake.output.metrics, sep='\t',
                   index=False, na_rep='NA')


if __name__ == '__main__':
    main()
