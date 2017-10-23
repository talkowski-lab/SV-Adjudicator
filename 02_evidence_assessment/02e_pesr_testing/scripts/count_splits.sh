#!/bin/bash
#
# split_stack.sh
#
# 
#
# Copyright (C) 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.


time_limit=7m
max_reads=1000000
bam_indexes="bam_indexes"
rerun_limit=7

function usage() {
cat <<EOF
Usage: count_splits.sh [-h] [-t TIME_LIMIT] [-r READ_LIMIT]
                       [-d BAM_INDEXES] [-n RERUN_LIMIT]
                       <sample_list> <count_table> <log> [region ...]

Query SR evidence at a breakpoint.

Required arguments:
  sample_list       List of samples to query
  count_table       Output table of split counts
  log               AWS query log
  region            Region(s) to query in each sample

Optional arguments:
  -t TIME_LIMIT     Time limit for split assembly. [4m]
  -r READ_LIMIT     Read count limit for PE check. [1000000]
  -d BAM_INDEXES    Directory of downloaded BAM indexes. [bam_indexes]
  -n RERUN_LIMIT    Number of times to attempt rerun upon AWS failure. [5]
EOF
}

if [[ $# -eq 0 ]]; then
  usage
  exit 0
fi

while getopts t:r:d:n:h opt; do
  case $opt in 
    t)
      time_limit=$OPTARG
      ;;
    r)
      max_reads=$OPTARG
      ;;
    d)
      bam_indexes=$OPTARG
      ;;
    n)
      rerun_limit=$OPTARG
      ;;
    h)
      usage
      exit 0
      ;;
  esac
done
shift $(( OPTIND - 1 ))

sample_list=$(readlink -f $1)
fout=$(readlink -f $2)
log=$(readlink -f $3)
regions="${@:4}"

# Use absolute path after moving to index subdirectory
count_splits=$(readlink -f ./scripts/count_splits.py)
bam_indexes=$(readlink -f $bam_indexes)

# Work from bam index directory to eliminate download time
cd bam_indexes

# Test each sample associated with variant
while read sample; do
  bam=$(s3bam $sample)

  # Check that bam path is valid
  aws s3 ls $bam &> /dev/null
  if [[ $? -ne 0 ]]; then
    echo "ERROR: $bam does not exist"
    exit 1
  fi

  # Track status
  aws_success=1
  rerun=0
  logged=0

  # Try to pull reads from S3 bucket until success or attempt limit reached
  while [[ $aws_success -ne 0 ]]; do
    # Give up after specified number of tries
    if [[ $rerun -eq "${rerun_limit}" ]]; then
      echo -e "RERUN\t${sample}" >> $log
      logged=1
      break
    fi

    # Write to stdout; write to file at end of loop
    samtools view -h $bam ${regions} 2> /dev/null \
      | timeout $time_limit $count_splits $sample 2> /dev/null

    # Check if S3 access succeeded
    pstat=( ${PIPESTATUS[*]} )
    aws_success=${pstat[0]}
    split_success=${pstat[1]}
 
    # Abort if read count exceeded
    if [[ $split_status -eq 2 ]]; then
      echo -e "READ_LIMIT\t${sample}" >> $log
      logged=1
      break
    fi
    
    # Abort if time limit exceeded
    if [[ $split_status -eq 124 ]]; then
      echo -e "TIMEOUT\t${sample}" >> $log
      logged=1
      break
    fi
 
    # Tick the attempt counter 
    rerun=$(($rerun+1))
  done
  
  if [[ $logged -eq 0 ]]; then
    echo -e "SUCCESS\t${sample}" >> $log
  fi
done < <(sed -e '/^name/d' $sample_list | cut -f2) > $fout
