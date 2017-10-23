#!/usr/bin/env python
import argparse
from collections import deque
import numpy as np
import pandas as pd
import pysam
import boto3
import sys
def load_s3vcf(bucket, vcf_path, index_filename=None):
    """
    Load an S3-hosted VCF.

    Parameters
    ----------
    bucket : str
        S3 bucket
    vcf_path : str
        S3 key

    Returns
    -------
    vcf : pysam.VariantFile
    """
    s3 = boto3.client('s3')
    url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={'Bucket': bucket, 'Key': vcf_path},
            ExpiresIn=86400)

    return pysam.VariantFile(url, index_filename=index_filename)
def filter_records(record):
    """
    Filter VCF records to those informative for BAF genotyping.

    Returns only records which match all of the following criteria:
    1) Biallelic
    2) SNP
    3) FILTER == PASS

    Parameters
    ----------
    records : iterator of pysam.VariantRecords

    Returns
    ------
    record : pysam.VariantRecord
    """

    # for record in records:
    # Restrict to biallelic sites
    if len(record.alleles) > 2:
        return

    # Restrict to variants which PASS
    if record.filter.keys() != ['PASS']:
        return

    # Restrict to SNPs
    ref, alt = record.alleles
    if len(ref) > 1 or len(alt) > 1:
        return

    return record
def calc_BAF(record, samples=None):
    """

    Parameters
    ----------
    record : pysam.VariantRecord
    samples : list of str, optional
        Subset of samples in record to consider

    Returns
    -------
    bafs : np.ndarray of np.float
        BAF at site for each sample
    """

    def _is_het(sample):
        return record.samples[sample]['GT'] == (0, 1)

    def _calc_BAF(sample):
        if not _is_het(sample):
            return np.nan

        DP = record.samples[sample]['DP']
        AD = record.samples[sample]['AD']

        if DP!=None and DP > 10: # SNP sites with >10 DP are included in BAF profile
            return AD[0] / DP
        else:
            return np.nan

        

    if samples is None:
        samples = record.samples.keys()

    bafs = np.atleast_2d(np.array([_calc_BAF(sample) for sample in samples], dtype=np.float))

    return bafs
def normalize_bafs(bafs, max_std=0.2):
    """
    Normalize BAFs and exclude outlying sites
    Normalize so per variant median BAF==0.5. Ignore sites with more than 0.2 standard deviation across samples. 
    
    Parameters
    ----------
    bafs : np.ndarray (n_sites x n_samples)
    max_std : float, optional
        Maximium standard deviation permitted at a site

    Returns
    -------
    normalized_bafs : np.ndarray
    """

    # Convert to n_samples x n_sites
    bafs = bafs.transpose()

    # Remove variants not informative in any sample (all NA)
    bafs = bafs.loc[:, ~bafs.isnull().all()]  # .copy()

    # Remove sites with excessive variance
    # Permit sites with a single sample (SD=NA)
    std = bafs.std()
    bafs = bafs.loc[:, ((std < max_std) | std.isnull())]

    # Center each site's median BAF at 0.5
    bafs = bafs - bafs.median()
    bafs = bafs + 0.5

    return bafs
    
#################
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument('bed', help='GATK VCF.')
    # parser.add_argument('vcf', help='GATK VCF.')
    # parser.add_argument('--tbi', help='Path to VCF index. (Required for tabix '
                        # 'querying when accessing an S3-hosted VCF.)')
    # parser.add_argument('-w', '--window', type=int, default=None,
                        # help='Window outside variant start and end to query '
                        # 'for SNPs. Defaults to CNV length if not specified.')
    # parser.add_argument('-s', '--samples', type=str, default=None,
                    # help='Samples')
    # parser.add_argument('-i', '--ID', type=str, default='test',help='Samples')
    # parser.add_argument('-t', '--type', type=str, default='test',help='Samples')
    parser.add_argument('-b', '--batch', default='batch.txt')
                    # help='Samples')
    args = parser.parse_args()
    # if args.vcf.startswith('s3://'):
        # vcf_path = args.vcf[5:]
        # bucket = vcf_path.split('/')[0]
        # vcf_path = '/'.join(vcf_path.split('/')[1:])

        # vcf = load_s3vcf(bucket, vcf_path, args.tbi)

    # else:
    vcf = pysam.VariantFile(sys.stdin)
    # While loop to iterate over all records, then break if reach the end
    idmap = '/data/talkowski/Samples/SFARI/lists/519families_idmapping'
    idmap = pd.read_table(idmap, names='sample SSC'.split())
    while True:
        try:
            record=next(vcf)
            record=filter_records(record)
            if record:
                site=[record.pos]
                site=np.array(site, dtype=np.int)
                samples = list(vcf.header.samples)
                baf=calc_BAF(record)
                # print(baf.shape)
                baf = pd.DataFrame(baf)
                baf.columns = samples
                baf = baf.set_index(site)
                baf = normalize_bafs(baf)
                baf.index.name = 'sample'
                baf = baf.reset_index()
                bf = pd.melt(baf, id_vars=['sample'], var_name='pos', value_name='baf')
                bf = bf.loc[~bf.baf.isnull()]
                called_bafs = bf
                called_bafs['chrom'] = record.chrom
                called_bafs['pos'] = called_bafs.pos.astype(int)
                # called_bafs['name'] = record.chrom + '-' + called_bafs.pos.astype(str)
                # called_bafs['genotype'] = 'AB'
                cols = 'chrom pos baf sample'.split()
                called_bafs = called_bafs[cols]
                if not called_bafs.empty:
                    print(called_bafs)
                    # called_bafs = called_bafs.rename(columns=dict(sample='SSC')) ## these two lines are replacing IDs, no need if no need to chcange ID
                    # called_bafs = pd.merge(called_bafs, idmap, on='SSC', how='left') ## these two lines are replacing IDs, no need if no need to chcange ID
                    called_bafs[cols].to_csv('mytest1.snp', index=False, mode='a',header=False, sep='\t')
        except StopIteration:
            break
    
    
    # splist=[]
    # with open(args.batch,'r') as f:
        # for line in f:
            # splist.append(line.rstrip())
    # with open(args.bed,'r') as f:
        # for line in f:
            # if line[0]!="#":
                # dat=line.rstrip().split('\t')
                # chrom=dat[0]
                # start=int(dat[1])
                # end=int(dat[2])
                # id=dat[3]
                # samples=dat[4]
                # samplelist=samples.split(',')
                # type=dat[5]
        # chrom, pos = args.region.split(':')
        # start, end = [int(x) for x in pos.split('-')]

                # het_counts, called_bafs = preprocess(chrom, start, end, vcf)

            # Temporary testing output
                # idmap = '/data/talkowski/Samples/SFARI/lists/519families_idmapping'
                # idmap = pd.read_table(idmap, names='sample SSC'.split())

                # het_counts = het_counts.rename(columns=dict(sample='SSC'))
                # het_counts = pd.merge(het_counts, idmap, on='SSC', how='left')
                # cols = 'before inside after sample'.split()

                # het_counts[cols].to_csv('mytest.het', index=False, header=False, sep='\t')
                # called_bafs = called_bafs.rename(columns=dict(sample='SSC'))
                # called_bafs = pd.merge(called_bafs, idmap, on='SSC', how='left')
                # cols = 'name chrom pos baf genotype sample'.split()
                
                # called_bafs[cols].to_csv('mytest.snp', index=False, header=False, sep='\t')
                ###### Running BAF testing
                # Del=DeletionTest(het_counts[het_counts['sample'].isin(splist)],samplelist,end-start)
                # KS=KS2sample(called_bafs[called_bafs['sample'].isin(splist)],samplelist)
                # ks,ksp=KS.test(samplelist)
                # mean,delp=Del.Ttest(samplelist)
                # statis=Del.stats(samplelist)
                # print(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+id+'\t'+samples+'\t'+type+'\t'+str(mean)+','+str(delp)+"\t"+str(ks)+','+str(ksp)+'\t'+statis)


if __name__ == '__main__':
    main()
