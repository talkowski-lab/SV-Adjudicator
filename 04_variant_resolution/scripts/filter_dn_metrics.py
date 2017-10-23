#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import pandas as pd


def filter_metrics(metrics, cutoffs, svclass, source):
    """
    Filter metrics table based on cutoffs

    Arguments
    ---------
    metrics : pd.DataFrame
    cutoffs : pd.DataFrame
        Indexed by svclass and source. Columns are metrics, values are cutoffs.
        NA if no cutoff applied.
    svclass : str
        Valid types: SV, CNV, DEL, DUP, INV, BND
        SV = DEL DUP INV BND
        CNV = DEL DUP
    source : str
        Valid types: pesr, depth
    """

    cutoffs = cutoffs.loc[(svclass, source)].squeeze().dropna().to_dict()

    # Check svtypes
    if svclass == 'SV':
        svtypes = 'DEL DUP INV BND'.split()
    elif svclass == 'CNV':
        svtypes = 'DEL DUP'.split()
    else:
        svtypes = [svclass]
    sv_match = (metrics.svtype.isin(svtypes))

    # Check source
    if source == 'depth':
        source_match = (metrics.sources == 'depth')
    else:
        source_match = (metrics.sources != 'depth')

    if svclass == 'SV':
        filter_pass = False
        for metric, cutoff in cutoffs.items():
            if metric == 'svsize':
                continue
            filter_pass = filter_pass | (metrics[metric] >= cutoff)
    else:
        filter_pass = True
        for metric, cutoff in cutoffs.items():
            filter_pass = filter_pass & (metrics[metric] >= cutoff)

    return (filter_pass & sv_match & source_match)


def filter_dn_metrics(metrics, cutoffs):
    # Filter depth-only

    depth_del_rd_pass = filter_metrics(metrics, cutoffs, 'DEL', 'depth')
    depth_dup_rd_pass = filter_metrics(metrics, cutoffs, 'DUP', 'depth')
    pesr_rd_pass = filter_metrics(metrics, cutoffs, 'CNV', 'pesr')
    pesr_pass = filter_metrics(metrics, cutoffs, 'SV', 'pesr')

    rd_pass = (depth_del_rd_pass | depth_dup_rd_pass | pesr_rd_pass)

    metrics.insert(5, 'pesr_pass', pesr_pass)
    metrics.insert(5, 'rd_pass', rd_pass)

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics')
    parser.add_argument('cutoffs')
    parser.add_argument('fout')
    args = parser.parse_args()

    metrics = pd.read_table(args.metrics)

    cutoffs = pd.read_table(args.cutoffs)
    cutoffs = cutoffs.drop('direction', axis=1)

    cutoffs = cutoffs.pivot_table(index='svtype source'.split(),
                                  values='cutoff', columns='metric')

    filtered = filter_dn_metrics(metrics, cutoffs)
    filtered.to_csv(args.fout, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()
