#!/bin/bash

before=$1
after=$2

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR="/tmp/"
fi

echo "TMPDIR is $TMPDIR"

if [[ $before == *.gz ]]; then
beforeCatCommand="zcat"
else
beforeCatCommand="/bin/cat"
fi

if [[ $after == *.gz ]]; then
afterCatCommand="zcat"
else
afterCatCommand="/bin/cat"
fi

uuid1=$(uuidgen)
uuid2=$(uuidgen)
uuid3=$(uuidgen)
uuid4=$(uuidgen)
uuid5=$(uuidgen)
uuid6=$(uuidgen)


$beforeCatCommand $before |skipcomments|awk '$3=="exon"' |sortgff| bedtools merge > $TMPDIR/$uuid1
$afterCatCommand $after |skipcomments|awk '$3=="exon"' |sortgff| bedtools merge > $TMPDIR/$uuid2
beforeNts=$(cat $TMPDIR/$uuid1| awk '{print $3-$2}' |sum.sh)
afterNts=$(cat $TMPDIR/$uuid2| awk '{print $3-$2}' |sum.sh)
let diff=$beforeNts-$afterNts || true # "true" serves to avoid "let" exiting with status > 0 when its output is = 0
echo "nts before: $beforeNts"
echo "nts after: $afterNts"
echo "diff nts: $diff"
if [ ! $diff -eq 0 ]; then echo "ERROR: Nucleotide coverage differ before/after merging";
echo -e "###### Features present before and absent after:\n\n";
bedtools intersect -v -f 1 -a $before -b $after 2> /dev/null
echo -e "\n\n"
exit 1;
fi

# all transcript_ids should be present before/after
$beforeCatCommand $before |skipcomments| extractGffAttributeValue.pl transcript_id | sort|uniq > $TMPDIR/$uuid3
$afterCatCommand $after |skipcomments|extractGffAttributeValue.pl contains |sed 's/,/\n/g' | sort|uniq > $TMPDIR/$uuid4
diffLines=$(diff -q $TMPDIR/$uuid3 $TMPDIR/$uuid4 |wc -l)
echo "difflines: $diffLines"
if [ ! $diffLines -eq 0 ]; then echo "ERROR: List of transcripts differ before/after"; exit 1; fi

#intron list should be the same before/after:
$beforeCatCommand $before |skipcomments|sortgff| makeIntrons.pl - |awk '{print $1"_"$4"_"$5"_"$7}' |sort|uniq > $TMPDIR/$uuid5
$afterCatCommand $after |skipcomments|sortgff| makeIntrons.pl - |awk '{{print $1"_"$4"_"$5"_"$7}}' |sort|uniq > $TMPDIR/$uuid6
diffIntronLines=$(diff -q $TMPDIR/$uuid5 $TMPDIR/$uuid6 |wc -l)

echo "diffintronlines: $diffIntronLines"
if [ ! $diffIntronLines -eq 0 ]; then echo "ERROR: List of introns differ before/after" exit 1; fi

echo XXdoneXX
rm -f $TMPDIR/$uuid1
rm -f $TMPDIR/$uuid2
rm -f $TMPDIR/$uuid3
rm -f $TMPDIR/$uuid4
rm -f $TMPDIR/$uuid5
rm -f $TMPDIR/$uuid6
