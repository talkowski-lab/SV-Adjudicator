#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import os
import boto3
import io


class S3File:
    def __init__(self, bucket, key, s3=None, chunk_size=4096):
        self.bucket = bucket
        self.key = key
        self.name = os.path.basename(key)

        if s3 is None:
            self.s3 = boto3.resource('s3')
        else:
            self.s3 = s3

        obj = self.s3.Object(self.bucket, self.key)
        body = obj.get()['Body']

        self.reader = iter(lambda: body.read(amt=chunk_size), b'')

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return next(self.reader)

    def fileno(self):
        """dummy for pysam"""
        return 1024

    @property
    def closed(self):
        return False

def iterable_to_stream(iterable, buffer_size=io.DEFAULT_BUFFER_SIZE):
    class IterStream(io.RawIOBase):
        def __init__(self):
            self.leftover = None
        def readable(self):
            return True
        def readinto(self, b):
            try:
                l = len(b)
                chunk = self.leftover or next(iterable)
                output, self.leftover = chunk[:l], chunk[l:]
                b[:len(output)] = output
                return len(output)
            except StopIteration:
                return 0
    return io.BufferedReader(IterStream(), buffer_size=buffer_size)


def s3iter(bucket, key, chunk_size=4096):
    s3 = boto3.resource('s3')
    obj = s3.Object(bucket, key)
    body = obj.get()['Body']

    for chunk in iter(lambda: body.read(amt=chunk_size), b''):
        yield chunk


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    args = parser.parse_args()

    s3file()


if __name__ == '__main__':
    main()
