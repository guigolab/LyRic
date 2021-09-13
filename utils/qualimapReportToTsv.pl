#!/usr/bin/env perl
use strict;
use warnings;

open F, "$ARGV[0]" or die $!;

my $mapped_bases=0;
my $insertions=0;
my $deletions=0;
my $mismatches=0;
my $bamfile;

while(<F>){
	chomp;
	if($_=~/\s+number of mapped bases = (\S+) bp/){
		$mapped_bases=$1;
		$mapped_bases=~s/,//g;
	}
	elsif($_=~/\s+bam file = (\S+)/){
		$bamfile=$1;
	}
	elsif($_=~/\s+number of mismatches = (\S+)/){
		$mismatches=$1;
		$mismatches=~s/,//g;
	}
	elsif($_=~/\s+number of insertions = (\S+)/){
		$insertions=$1;
		$insertions=~s/,//g;
	}
	elsif($_=~/\s+number of deletions = (\S+)/){
		$deletions=$1;
		$deletions=~s/,//g;
	}
}

my $insertionsPerMappedBase=$insertions/$mapped_bases;
my $deletionsPerMappedBase=$deletions/$mapped_bases;
my $mismatchesPerMappedBase=$mismatches/$mapped_bases;
my $globalErrorRate=($insertions+$deletions+$mismatches)/$mapped_bases;

print "$bamfile\tmismatchesPerMappedBase\t$mismatchesPerMappedBase
$bamfile\tdeletionsPerMappedBase\t$deletionsPerMappedBase
$bamfile\tinsertionsPerMappedBase\t$insertionsPerMappedBase
$bamfile\tglobalErrorRate\t$globalErrorRate
";
