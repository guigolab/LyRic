#!/bin/bash

set -e

source ~/julien_utils/check_pipe_exit.sh

fileList=$1
#bigBedAutoSqlFile=$2
chrFile=$2
outDir=$3
cat $fileList| while read path attrs; do
bn=`basename $path .bed`
echo -e "####### $path\n\n"
#check if bed6 or bed12 format
bedFields=`cat $path | head -n1 | wc -w`
if [ "$bedFields" == "12" ]; then 
bigBedAutoSqlFile="/users/rg/jlagarde/bed12.as"
else
if [ "$bedFields" == "6" ]; then 
bigBedAutoSqlFile="/users/rg/jlagarde/bed6.as"
else
if [ "$bedFields" == "9" ]; then 
bigBedAutoSqlFile="/users/rg/jlagarde/bed9.as"
else 
echoerr "WRONG FORMAT";
exit 1;
fi
fi
fi

new="$outDir/$bn.bb"
mkdir -p $outDir

if [ ! -f $new ];
then
echo "Sorting $path..." >&2
jcat $path | sort -k1,1 -k2,2n > /tmp/$bn.sort
if [ "$(check_exit)" != "0" ]; then echo "ERROR" >&2; exit 1; fi
echo "done" >&2
bedToBigBed -as=$bigBedAutoSqlFile -type=bed$bedFields -extraIndex=name /tmp/$bn.sort $chrFile $new
echo "done. Output in $new ." >&2
rm /tmp/$bn.sort
else
echo "$new already exists. No need to re-generate it" >&2
fi
printf "$new\t$attrs\n"
done
