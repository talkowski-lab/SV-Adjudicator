#!bin/bash
# This script will merge split rdtest statistics by split number
# The output will be there be 1 file for each chrom+batch+source for allosomes
# For sex chromosome the result will be one file per sex chromosome for each sex
# The usage format is ```bash rdtest_mergesplit.sh {batch} {source} {chrom}
batch=$1
source=$2
chrom=$3
inputfolder="split_srtest"
outputfolder="srtest"
input=()
for item in  $inputfolder/$batch.$source.$chrom.*.stats; do
        input+=($item)
done

if [ "$chrom" = "X" ]; then
    cat $inputfolder/$batch.$source.X.*.females.stats \
          | sed -r -e '/^name/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $inputfolder/$batch.$source.X.females.stats
    cat $inputfolder/$batch.$source.X.*.males.stats \
          | sed -r -e '/^name/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $inputfolder/$batch.$source.X.males.stats
    python scripts/sr_merge_allosomes.py $batch $source X $inputfolder $outputfolder
elif [ "$chrom" = "Y" ]; then
    cat $inputfolder/$batch.$source.Y.*.females.stats \
          | sed -r -e '/^name/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $inputfolder/$batch.$source.Y.females.stats
    cat $inputfolder/$batch.$source.Y.*.males.stats \
          | sed -r -e '/^name/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $inputfolder/$batch.$source.Y.males.stats
    python scripts/sr_merge_allosomes.py $batch $source Y $inputfolder $outputfolder
else
    cat $inputfolder/$batch.$source.$chrom.*.stats \
          | sed -r -e '/^name/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $outputfolder/$batch.$source.$chrom.stats
fi

