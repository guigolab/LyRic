#!/usr/bin/env perl

use strict;
use warnings;

my $simpleGffCompare=$ARGV[0];
my $bed=$ARGV[1];


open TSV, $simpleGffCompare or die "$!";

my %txToClass=();
while(<TSV>){
	chomp;
	my @line=split "\t";
	die "Wrong format for file $simpleGffCompare.\n" unless ($#line==1);
	$txToClass{$line[0]}=$line[1];
}
close TSV;

open BED, $bed or die "$!";

while(<BED>){
	chomp;
	my @line=split "\t";
	if (exists $txToClass{$line[3]}){
		if($txToClass{$line[3]} eq 'Equal' || $txToClass{$line[3]} eq 'Included'){
			$line[8]="0,0,0";
		}
		elsif($txToClass{$line[3]} eq 'Extends' || $txToClass{$line[3]} eq 'Overlaps' || $txToClass{$line[3]} eq 'Antisense' || $txToClass{$line[3]} eq 'Intronic' || $txToClass{$line[3]} eq 'Intergenic'){
			$line[8]="180,130,0";

		}
		else{
			die "Unknown category for transcript $line[3]: $txToClass{$line[3]}.\n";
		}

	}
	else{
		die "$line[3] absent from $simpleGffCompare.\n";
	}
	print join("\t", @line)."\n";
}
