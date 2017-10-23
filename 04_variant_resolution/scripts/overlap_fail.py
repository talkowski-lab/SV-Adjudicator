#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import pysam
from svtools.vcfcluster import VCFCluster


def is_batch_only(cluster, batch):
    ids = [s.record.id for s in cluster.records]
    return all([batch in x for x in ids])


def get_sources(header):
    for record in header.records:
        rec = str(record)
        if rec.startswith('##source='):
            return rec.strip().split('=')[1].split(',')

    return []


def overlap_fail(phase1, pilot, header, dist=300, frac=0.1):
    svc = VCFCluster([phase1, pilot], dist=dist, frac=frac, preserve_ids=True)
    sources = get_sources(header)

    for cluster in svc.cluster(merge=False):
        # skip phase1-only variants and pilot variants that overlap with
        # rejected phase1 variants
        if is_batch_only(cluster, 'Pilot'):
            # Make new record and merge in pilot samples
            record = header.new_record()
            record = cluster.merge_record_data(record)
            record = cluster.merge_record_formats(record, sources,
                                                  call_sources=True)
            record.info['MEMBERS'] = [r.record.id for r in cluster.records]

            # check that we're not overclustering Pilot variants
            if len(cluster.records) != 1:
                raise Exception('Multiple Pilot variants clustered')

            record.id = cluster.records[0].record.id
            yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phase1', help='Failed Phase1 calls')
    parser.add_argument('pilot', help='Filtered pilot calls')
    parser.add_argument('fout')
    args = parser.parse_args()

    phase1 = pysam.VariantFile(args.phase1)
    pilot = pysam.VariantFile(args.pilot)

    header = phase1.header.copy()
    for sample in pilot.header.samples:
        header.add_sample(sample)

    # TEMPORARY; not all vcfs were clustered since adding MEMBERS
    if 'MEMBERS' not in header.info.keys():
        info = ('##INFO=<ID=MEMBERS,Number=.,Type=String,'
                'Description="IDs of cluster\'s constituent records.">')
        header.add_line(info)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')
    fout = pysam.VariantFile(fout, mode='w', header=header)

    for record in overlap_fail(phase1, pilot, header):
        fout.write(record)


if __name__ == '__main__':
    main()
