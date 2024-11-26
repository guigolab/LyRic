#!/bin/bash

shopt -s expand_aliases
source ~/.bash_profile
cat $1 | \
#centering coordinates:
perl -F"\t" -lane 'unless ( $F[1] =~ /^\d+$/ && $F[2]=~ /^\d+$/ ){die "Input is not in BED format. Cannot continue.\n"} if ( ( $F[2] - $F[1] ) % 2 == 1 ) {
            $F[2]++;
        } $midpoint = ( $F[1] + $F[2] ) / 2;; $F[1]= $midpoint -1; $F[2] = $midpoint; $F[3].=".$midpoint"; print join ("\t", @F)' | \
#extending +/- $2 nts on each side using bedtools
bedtools slop -b $3 -i stdin -g $2  | \
#calculating coverage of each base
bedtools coverage -s -d -a stdin -b $4 | \
# adjust coordinates to -$2, +$2
perl -F"\t" -slane '$F[3]=~/.+\.(\d+)$/; $midpoint=$1; $F[3]=~s/(.+)\.\d+$/$1/; $centeredStart=($F[1]-$midpoint)+$F[6]; $centeredStart=-$centeredStart if($F[5] eq "-"); $F[6]=$centeredStart; print "$F[3]\t$F[6]\t$F[7]"'

