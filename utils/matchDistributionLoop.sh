#!/bin/bash

pThreshold=0.05
passes=$1 # int
doKStest=$2 # 1/0 boolean, do KS test or not
targetDist=$3
subjectDist=$4
bins=$5
breakIfKSTest=$6 # 1/0 boolean, break loop if KS test gives p>$pThreshold

ln -f -s $subjectDist $subjectDist.sampled.0.tsv
for i in `seq 1 $passes`; do
let j=$i-1
echo "Current pass: $i"
matchDistribution.pl $targetDist $bins $subjectDist.sampled.$j.tsv > $subjectDist.sampled.$i.tsv
if [ $doKStest == 1 ]; then 
outR="$(matchDistributionKStest.r $targetDist $subjectDist.sampled.$i.tsv)"
ksTestBool=`echo $outR| awk -v p=$pThreshold '{if($1>p) print 1; else print 0}'`
echo "p-value (Kolmogorov-Smirnov test): $outR"
else
ksTestBool=0
fi
rm -f $subjectDist.sampled.$j.tsv
if [ $breakIfKSTest == 1 ]; then
if [ $ksTestBool == 1 ]; then 
echo "p-value > 0.05, i.e. the two distribution now look reasonably similar. Quitting iterations after $i passes. Final result is in $subjectDist.sampled.$i.tsv"
echo
break;
fi
fi
done 
