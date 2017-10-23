
configfile: 'config.yaml'

subworkflow background:
    workdir: "."
    snakefile: "background.snake"

from collections import deque
import os
import numpy as np
import pandas as pd

with open(config['bed']) as bedfile:
    next(bedfile)
    NAMES = [line.strip().split()[3] for line in bedfile]
#NAMES = ['polymorphic_cnv_1858188']

PREFIX = os.path.splitext(os.path.basename(config['bed']))[0]
REGIONS = {}
with open(background('calls/{prefix}.sr_windows.txt'.format(prefix=PREFIX))) as regionfile:
    for line in regionfile:
        name, svtype, regions = line.strip().split('\t')
        if name in NAMES:
            REGIONS[name] = regions


rule all:
    input:
        'data/{prefix}.merged_pvals.txt'.format(prefix=PREFIX),
        'data/{prefix}.variant_stats.txt'.format(prefix=PREFIX)
        
rule count_splits:
    input:
        sample_list=background('sample_lists/{name}.txt')
    output:
        counts='split_counts/{name}.txt',
        log='count_logs/{name}.log'
    params:
        regions=lambda wildcards: REGIONS[wildcards.name],
    shell:
        """
        ./count_splits.sh {input.sample_list} {output.counts} {output.log} "{params.regions}"
        """

rule collapse_splits:
    input:
        bed=config['bed'],
        samples='sample_lists/{name}.txt',
        counts='split_counts/{name}.txt'
    params:
        window=config['window']
    output:
        'coord_counts/{name}.txt'
    script:
        "combine_splits.py"

rule calc_pvals:
    input:
        'coord_counts/{name}.txt'
    output:
        'pvals/{name}.txt'
    script:
        "srtest.py"

rule get_variant_stats:
    input:
        config['bed']
    output:
        'data/{prefix}.variant_stats.txt'
    run:
        bed = pd.read_table(input[0])
        bed['svsize'] = bed.end - bed.start
        bed['log_svsize'] = np.log10(bed.svsize)
        bed['num_samples'] = bed.samples.str.split(',').apply(len)
        bed['vaf'] = bed.num_samples / 2076
        cols = 'name svtype batch log_svsize vaf'.split()
        bed[cols].to_csv(output[0], index=False, sep='\t')
        
rule combine_pvals:
    input:
        pvals=expand('pvals/{name}.txt', name=NAMES),
        stats='data/{prefix}.variant_stats.txt'
    output:
        'data/{prefix}.merged_pvals.txt'
    run:
        pvals = deque()
        for pval in input.pvals:
            pval = pd.read_table(pval)
            pvals.append(pval)
        pvals = pd.concat(pvals)
        pvals.to_csv(output[0], sep='\t', index=False)


