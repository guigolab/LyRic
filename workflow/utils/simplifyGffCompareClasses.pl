#!/usr/bin/env perl

use strict;
use warnings;

open F, "$ARGV[0]" or die $!;

while(<F>){
	my @line=split "\t";
	my $gffCompareClass=$line[3];
	my $subject=$line[4];
	my @subj=split(" ", $subject);
	die "Can't process multi-sample gffcompare output.\n" if ($#subj>0);
	$subject=~/q1:(\S+?)\|(\S+?)\|\S+/;
	$subject=$2;
	if($gffCompareClass eq '='){
		$gffCompareClass="Equal";
	}
	elsif($gffCompareClass eq 'c'){
		$gffCompareClass="Included"
	}
	elsif($gffCompareClass eq 'k'){
		$gffCompareClass="Extends"
	}
	elsif($gffCompareClass eq 'j'){
		$gffCompareClass="Overlaps"
	}
	elsif($gffCompareClass eq 'e'){
		$gffCompareClass="Overlaps"
	}
	elsif($gffCompareClass eq 'o'){
		$gffCompareClass="Overlaps"
	}
	elsif($gffCompareClass eq 'm'){
		$gffCompareClass="Overlaps"
	}
	elsif($gffCompareClass eq 'n'){
		$gffCompareClass="Overlaps"
	}
	elsif($gffCompareClass eq 's'){
		$gffCompareClass="Antisense"
	}
	elsif($gffCompareClass eq 'x'){
		$gffCompareClass="Antisense"
	}
	elsif($gffCompareClass eq 'i'){
		$gffCompareClass="Intronic"
	}
	elsif($gffCompareClass eq 'y'){
		$gffCompareClass="Intergenic"
	}
	elsif($gffCompareClass eq 'p'){
		$gffCompareClass="Intergenic"
	}
	elsif($gffCompareClass eq 'u'){
		$gffCompareClass="Intergenic"
	}
	else{
		die "Unknown class '$gffCompareClass'\n"
	}
	$line[3]=$gffCompareClass;
	print "$subject\t$gffCompareClass\n"

}
