#!/usr/bin/env python
# -*- coding: utf-8 -*-
#freq
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
from collections import deque, defaultdict
import numpy as np
import scipy.stats as ss
import pandas as pd
import pysam
import svtools.utils as svu


def process_rdtest(rdtest):
    """Standardize rdtest column names"""

    # Drop metadata columns (available from VCF) and rename CNVID
    skip_cols = 'chr Start End SampleIDs Type'.split()
    rdtest = rdtest.drop(skip_cols, axis=1).rename(columns={'CNVID': 'name'})

    numeric_cols = 'Median_Power P 2ndMaxP Median_Rank Median_Separation'
    numeric_cols = numeric_cols.split()

    # Replace strings with NA
    for col in numeric_cols:
        repl = ['All_samples_called_CNV_no_analysis',
                'No_samples_for_analysis']
        rdtest[col] = rdtest[col].replace(repl, np.nan).astype(np.float)

    rdtest['log_pval'] = -np.log10(rdtest.P)
    rdtest['log_2ndMaxP'] = -np.log10(rdtest['2ndMaxP'])

    maxp = rdtest.loc[rdtest.log_pval != np.inf, 'log_pval'].max()
    max2p = rdtest.loc[rdtest.log_2ndMaxP != np.inf, 'log_2ndMaxP'].max()

    rdtest.loc[rdtest.log_pval == np.inf, 'log_pval'] = maxp + 5
    rdtest.loc[rdtest.log_2ndMaxP == np.inf, 'log_2ndMaxP'] = max2p + 5

    rdtest.log_pval = rdtest.log_pval.abs()
    rdtest.log_2ndMaxP = rdtest.log_2ndMaxP.abs()

    return rdtest


def process_srtest(srtest):
    metrics = 'log_pval called_median bg_median'.split()

    # remove -0.0 (temporary, should fix in SR-test)
    srtest.log_pval = srtest.log_pval.abs()

    # force one-sided (temporary, should fix in SR-test)
    srtest.loc[srtest.bg_median > srtest.called_median, 'log_pval'] = 0

    srtest = srtest.pivot_table(index='name', values=metrics, columns='coord')
    srtest.columns = ['_'.join(col[::-1]).strip()
                      for col in srtest.columns.values]
    srtest = srtest.reset_index()

    return srtest


def process_petest(petest):
    # remove -0.0 (temporary, should fix in PE-test)
    petest.log_pval = petest.log_pval.abs()

    # force one-sided (temporary, should fix in PE-test)
    petest.loc[petest.bg_median > petest.called_median, 'log_pval'] = 0

    return petest


def process_baftest(baftest):
    skip_cols = 'chrom start end samples svtype'.split()
    baftest = baftest.drop(skip_cols, axis=1)

    baftest['KS_log_pval'] = (- np.log10(baftest.KS_pval)).abs()
    baftest['del_loglik'] = -baftest.del_loglik

    repl = 'Potential ROHregion or reference error'
    baftest.delstat = baftest.delstat.replace(repl, 'Ref_error')

    return baftest


def preprocess(df, dtype):
    if dtype == 'RD':
        return process_rdtest(df)
    elif dtype == 'SR':
        return process_srtest(df)
    elif dtype == 'PE':
        return process_petest(df)
    elif dtype == 'BAF':
        return process_baftest(df)
    else:
        return df


def _is_parent(s):
    return s.endswith('fa') or s.endswith('mo')


def _is_child(s):
    return s.endswith('p1') or s.endswith('s1')


def get_inh_rate(called):
    quads = defaultdict(list)

    for sample in called:
        quad, member = sample.split('.')
        quads[quad].append(member)

    n_called = len([s for s in called if _is_child(s)])

    n_inh = 0
    for quad, members in quads.items():
        if 'fa' in members or 'mo' in members:
            if 'p1' in members:
                n_inh += 1
            if 's1' in members:
                n_inh += 1

    return n_inh / n_called


def process_metadata(variants, bed=False, batch_list=None):
    if bed:
        samples = [s.strip() for s in batch_list.readlines()]
    else:
        samples = list(variants.header.samples)

    parents = [s for s in samples if _is_parent(s)]
    children = [s for s in samples if _is_child(s)]
    n_parents = len(parents)
    n_children = len(children)

    metadata = deque()
    for variant in variants:
        # bed record
        if bed:
            if variant.startswith('#'):
                continue
            data = variant.strip().split()
            called = data[4].split(',')
            name = data[3]
            svtype = data[5]
        # VCF record
        else:
            called = svu.get_called_samples(variant)
            name = variant.id
            svtype = variant.info['SVTYPE']

        # Calculate parental VF
        parents = [s for s in called if _is_parent(s)]
        parental_vf = len(parents) / n_parents

        children = [s for s in called if _is_child(s)]
        child_vf = len(children) / n_children

        if child_vf > 0:
            inh_rate = get_inh_rate(called)
        else:
            inh_rate = 0

        dat = [name, svtype, parental_vf, child_vf, inh_rate]
        metadata.append(dat)

    metadata = np.array(metadata)
    cols = 'name svtype parental_vf child_vf inh_rate'.split()
    metadata = pd.DataFrame(metadata, columns=cols)
    return metadata


def add_pesr(evidence):
    evidence['PESR_called_median'] = (evidence['PE_called_median'] +
                                      evidence['SR_sum_called_median'])
    evidence['PESR_bg_median'] = (evidence['PE_bg_median'] +
                                  evidence['SR_sum_bg_median'])

    def calc_p(row):
        pval = ss.poisson.cdf(row.PESR_bg_median, row.PESR_called_median)
        return np.abs(-np.log10(pval))

    evidence['PESR_log_pval'] = evidence.apply(calc_p, axis=1)

    one_sided_mask = (evidence.PESR_bg_median > evidence.PESR_called_median)
    evidence.loc[one_sided_mask, 'PESR_log_pval'] = 0

    return evidence


def make_columns():
    PE_names = ('log_pval called_median bg_median').split()
    PESR_names = ['PESR_' + name for name in PE_names]
    PE_names = ['PE_' + name for name in PE_names]

    SR_names = ('posA_log_pval posB_log_pval sum_log_pval '
                'posA_called_median posB_called_median sum_called_median '
                'posA_bg_median posB_bg_median sum_bg_median').split()
    SR_names = ['SR_' + name for name in SR_names]

    BAF_names = ('delstat snp_ratio del_loglik dupstat KS_stat KS_log_pval '
                 'total_case_snps total_snps n_nonROH_cases n_samples '
                 'mean_control_snps n_nonROH_controls n_controls').split()
    BAF_names = ['BAF_' + name for name in BAF_names]

    RD_names = ('Median_Power P 2ndMaxP Model Median_Rank Median_Separation '
                'log_pval log_2ndMaxP').split()
    RD_names = ['RD_' + name for name in RD_names]

    metadata_names = 'name svtype parental_vf child_vf inh_rate'.split()

    return (metadata_names + PE_names + SR_names + PESR_names + RD_names +
            BAF_names)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--variants', required=True, help='Default VCF')
    parser.add_argument('-r', '--RDtest')
    parser.add_argument('-b', '--BAFtest')
    parser.add_argument('-s', '--SRtest')
    parser.add_argument('-p', '--PEtest')
    parser.add_argument('--batch-list', type=argparse.FileType('r'))
    parser.add_argument('-d', '--bed', action='store_true', default=False)
    parser.add_argument('fout')
    args = parser.parse_args()

    if args.bed:
        if not hasattr(args, 'batch_list'):
            raise Exception('batch list must be specified when passing a bed')
        variants = open(args.variants)
        dtypes = 'RD BAF'.split()
    else:
        variants = pysam.VariantFile(args.variants)
        dtypes = 'PE SR RD BAF'.split()

    metadata = process_metadata(variants, args.bed, args.batch_list)
    metadata = metadata.set_index('name')

    evidence = deque()

    BAF_names = ('chrom start end name samples svtype delstat snp_ratio '
                 'del_loglik dupstat KS_stat KS_pval total_case_snps '
                 'total_snps n_nonROH_cases n_samples mean_control_snps '
                 'n_nonROH_controls n_controls').split()

    for dtype in dtypes:
        dtable = getattr(args, dtype + 'test')
        if dtable is None:
            continue

        names = BAF_names if dtype == 'BAF' else None
        df = pd.read_table(dtable, names=names)

        df = preprocess(df, dtype)
        df = df.rename(columns=lambda c: dtype + '_' + c if c != 'name' else c)
        df = df.set_index('name')
        evidence.append(df)

    evidence = list(evidence)
    evidence = metadata.join(evidence, how='outer', sort=True)
    evidence = evidence.reset_index().rename(columns={'index': 'name'})

    has_petest = (getattr(args, 'PEtest') is not None)
    has_srtest = (getattr(args, 'SRtest') is not None)
    if not args.bed and has_petest and has_srtest:
        evidence = add_pesr(evidence)

    # Replace infinite log-pvals
    LOG_CEIL = 300
    evidence = evidence.replace(np.inf, LOG_CEIL)

    evidence = evidence.reindex(columns=make_columns())
    evidence.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
