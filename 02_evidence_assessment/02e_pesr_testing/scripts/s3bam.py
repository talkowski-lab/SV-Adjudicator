#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 ec2-user <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""

"""

import os
import boto3
import pysam

def load_s3bam(bucket, bam_path, filepath_index=None):
    s3 = boto3.client('s3')
    url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={'Bucket': bucket, 'Key': bam_path},
            ExpiresIn=86400)

    return pysam.AlignmentFile(url, filepath_index=filepath_index)


def load_bam(bam_path, index_dir=None):
    if bam_path.startswith('s3://'):
        s3path = bam_path[5:]
        bucket = s3path.split('/')[0]
        bam_path = '/'.join(s3path.split('/')[1:])

        # Get index if possible
        bam_name = os.path.basename(bam_path)
        if index_dir is None:
            filepath_index = None
        else:
            idx1 = bam_name + '.bai'
            idx2 = os.path.splitext(bam_name)[0] + '.bai'
            if os.path.exists(idx1):
                filepath_index = idx1
            elif os.path.exists(idx2):
                filepath_index = idx2
            else:
                filepath_index = None

        bam = load_s3bam(bucket, bam_path, filepath_index)
    else:
        bam = pysam.AlignmentFile(bam_path)

    return bam
