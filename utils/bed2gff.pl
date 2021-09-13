#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
#takes only BED6 as input, one exon per line

#awk '{print $1"\t.\t.\t"$2+1"\t"$3"\t.\t.\t.\t"$1"_"$2+1"_"$3}' "$@"

while (<STDIN>){
	if ($_=~/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\n/){
		my $chr=$1;
		my $start=$2+1;
		my $end=$3;
		my $transcript_id=$4;
		my $score=$5;
		my $strand=$6;
		print "$chr\tsource\texon\t$start\t$end\t$score\t$strand\t.\ttranscript_id \"$transcript_id\";\n";
	}
	else{
		die "Non BED6 line at line $. . Cannot continue: \n$_\n"
	}
}
