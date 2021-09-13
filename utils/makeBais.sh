#!/bin/bash

set -e

source ~/julien_utils/check_pipe_exit.sh

fileList=$1

for file in `cut -f1 $fileList`; do

echo "$file
" >&2
new=`dirname $file`/`basename $file`.bai

if [ ! -f $new ];
then
echo "Making BAI..." >&2
time=`now`
samtools index $file tmp.$time 
\mv tmp.$time $file.bai
echo "done. Output in $new ." >&2
else
echo "$new already exists. No need to re-generate it" >&2
fi
done
