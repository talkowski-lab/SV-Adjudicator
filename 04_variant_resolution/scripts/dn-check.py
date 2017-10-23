#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


import argparse
from collections import namedtuple, deque
import numpy as np
import pysam
import svtools.utils as svu
from svtools.cli.sr_test import SRBreakpoint
from svtools.cli.pe_test import PEBreakpoint


def check_record(record, fam_map, config, max_vaf=0.01):
    # get total number of parents for VF filtering
    n_parents = np.unique(fam_map[:, 1:]).shape[0]
    max_parents = max_vaf * n_parents

    # Convert genotypes to 0/1/2
    genotypes = [sum(filter(None, s['GT'])) for s in record.samples.values()]
    genotypes = np.array(genotypes, dtype=np.int64)

    # Pivot genotypes into (child, mother, father) tuples
    gt = np.take(genotypes, fam_map)

    # Restrict to rare variants
    called_parent = (gt[:, 1:] != 0)
    # parents appear once for each child
    n_called_parents = np.sum(called_parent.flatten()) / 2
    if n_called_parents >= max_parents:
        return

    # Restrict to variants called in children, then filter to calls with
    # no call in parent
    called = (gt[:, 0] != 0)
    no_parents = (gt[:, 1:].sum(axis=1) == 0)
    denovo = fam_map[np.where(called & no_parents)]

    if denovo.shape[0] != 0:
        samples = np.array(record.samples.keys(), dtype=str)
        # parents = np.unique(np.take(samples, denovo[:, 1:]))
        parents = np.unique(np.take(samples, denovo))
        dn_test(record, list(parents), config)


def dn_test(record, parents, config):
    called = svu.get_called_samples(record)

    for i, parent in enumerate(parents):
        others = parents[:i] + parents[i+1:]
        blacklist = called + others
        samples = record.samples.keys()
        whitelist = [s for s in samples if s not in blacklist]

        # PE Test
        pe = PEBreakpoint.from_vcf(record)
        pe.samples = [parent]

        pe.pe_test(whitelist, config.discfile,
                   n_background=160, window_in=50, window_out=500)
        stats = pe.stats
        stats['name'] = pe.name
        stats['sample'] = parent
        cols = 'name sample log_pval called_median bg_median'.split()
        stats[cols].to_csv(config.petest, sep='\t', index=False, header=False,
                           na_rep='NA')

        # SR Test
        sr = SRBreakpoint.from_vcf(record)
        sr.samples = [parent]

        sr.sr_test(whitelist, config.countfile,
                   n_background=160, window=50)

        pvals = sr.best_pvals
        pvals['sample'] = parent
        cols = 'name sample coord pos log_pval called_median bg_median'.split()
        pvals = pvals[cols].fillna(0)

        int_cols = ['pos']  # called_median bg_median'.split()
        for col in int_cols:
            pvals[col] = pvals[col].round().astype(int)
        pvals.log_pval = np.abs(pvals.log_pval)

        pvals.to_csv(config.srtest, sep='\t', index=False, header=False,
                     na_rep='NA')


class DnTestConfig:
    def __init__(self, petest, srtest, discfile, countfile):
        self.petest = petest
        self.srtest = srtest
        self.discfile = discfile
        self.countfile = countfile


def parse_fam(famfile, vcf):
    """
    Create mapping of children to parent indices in VCF header

    Returns
    -------
    fam_map : np.array([int, int, int])
        Indices of sample, mother, father in VCF header
    """

    samples = list(vcf.header.samples)
    idx_map = deque()

    names = 'fam sample father mother sex phenotype'.split()
    Sample = namedtuple('Sample', names)

    for line in famfile:
        data = line.strip().split()
        sample = Sample(*data)

        if sample.sample not in samples:
            continue
        # Require both parents to determine de novo status
        if sample.father != '0' and sample.mother != '0':
            indices = np.zeros(3)
            indices[0] = samples.index(sample.sample)
            indices[1] = samples.index(sample.mother)
            indices[2] = samples.index(sample.father)
            idx_map.append(indices)

    return np.array(idx_map, dtype=np.int64)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fam', type=argparse.FileType('r'))
    parser.add_argument('countfile')
    parser.add_argument('discfile')
    parser.add_argument('petest', type=argparse.FileType('w'), help='fout')
    parser.add_argument('srtest', type=argparse.FileType('w'), help='fout')
    args = parser.parse_args()

    countfile = pysam.TabixFile(args.countfile)
    discfile = pysam.TabixFile(args.discfile)

    config = DnTestConfig(args.petest, args.srtest, discfile, countfile)

    vcf = pysam.VariantFile(args.vcf)
    fam_map = parse_fam(args.fam, vcf)

    header = 'name sample log_pval called_median bg_median'.split()
    args.petest.write('\t'.join(header) + '\n')

    header = 'name sample coord pos log_pval called_median bg_median'.split()
    args.srtest.write('\t'.join(header) + '\n')

    for record in vcf:
        check_record(record, fam_map, config)


if __name__ == '__main__':
    main()
