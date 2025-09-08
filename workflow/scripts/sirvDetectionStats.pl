#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Data::Dumper;

my $sirvInfo=$ARGV[0];
my $gffCompareRefmap=$ARGV[1];

open F, "$sirvInfo" or die $!;

my %sirvConcentration=();
my %sirvLength=();

while (<F>){
	chomp;
	my @line=split "\t";
	$sirvConcentration{$line[0]}=$line[2];
	$sirvLength{$line[0]}=$line[1];
}


close F;

open G, "$gffCompareRefmap" or die $!;
my %sirvDetection=();
while (<G>){
	chomp;
	next if $_=~/^ref_gene_id\tref_id\tclass_code\tqry_id_list$/;
	my @line=split "\t";
	if($line[2] eq '='){
		$sirvDetection{$line[1]}="end-to-end"
	}
	elsif($line[2] eq 'c'){
		unless (exists $sirvDetection{$line[1]}){
			$sirvDetection{$line[1]}="partial"
		}
	}
	else{
		die "Wrong gffcompare refmap format: \n$_\n";
	}
}

foreach my $sirv (sort keys %sirvLength){
	my $detection;
	if(exists $sirvDetection{$sirv}){
		$detection=$sirvDetection{$sirv}
	}
	else{
		$detection="absent";
	}
	print "$sirv\t$sirvLength{$sirv}\t$sirvConcentration{$sirv}\t$detection\n";

}
