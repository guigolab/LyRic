#!/bin/bash

shopt -s expand_aliases
source ~/.bash_profile
cat $1 | \
#centering coordinates:
#perl -F"\t" -lane 'unless ( $F[1] =~ /^\d+$/ && $F[2]=~ /^\d+$/ ){die "Input is not in BED format. Cannot continue.\n"} if ( ( $F[2] - $F[1] ) % 2 == 1 ) {
#            $F[2]++;
#        } $midpoint = ( $F[1] + $F[2] ) / 2;; $F[1]= $midpoint -1; $F[2] = $midpoint; $F[3].=".$midpoint"; print join ("\t", @F)' | \
#extending +/- $2 nts on each side using bedtools
#bedtools slop -b $3 -i stdin -g $2  | \
#calculating coverage of each base
bedtools coverage -s -d -a stdin -b $3 |\
# adjust coordinates
perl -F"\t" -slane '$length=$F[2]-$F[1]; if($F[5] eq "+"){$F[6]=$F[6] / $length} elsif($F[5] eq "-"){$F[6]=($length - $F[6]) / $length} else{die} $F[6]=sprintf("%.2f", $F[6]); print "$F[3]\t$F[6]\t$F[7]"'
