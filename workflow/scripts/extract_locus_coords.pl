#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
#use Bio::DB::GFF;
#use Bio::SeqIO;
#use Data::Dumper;

#input: GFF/GTF
#output: BED

open I, "$ARGV[0]" or die $!;


my %transcripts_exonstarts;
my %transcripts_exonends;
my %transcripts_strand;
my %transcripts_chr;

while (<I>){
	if ($_=~/^(\S+)\t(\S+)\texon\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*gene_id "(\S+)?";.*\n/){
		my $chr=$1;
		my $source=$2;
		my $start=$3;
		my $end=$4;
		my $strand=$5;
		my $locus=$7;
		push(@{$transcripts_exonstarts{$locus}},$start);
		push(@{$transcripts_exonends{$locus}},$end);
		if((exists $transcripts_strand{$locus}) && ($strand ne $strand) ){
			print STDERR "WARNING: Locus $locus has transcripts on different strands!!! Locus $locus marked with strand '.' in output.\n";
			$strand='.';
		}
		$transcripts_strand{$locus}=$strand;
		$transcripts_chr{$locus}=$chr;

	}
}

foreach my $l (keys %transcripts_exonstarts){

		unless($#{$transcripts_exonstarts{$l}}==$#{$transcripts_exonends{$l}}){
		print STDERR "Nr of exon starts differs from nr of exon ends for locus $l!!!\n";
	}
	my $start_l='';
	my $end_l='';

	@{$transcripts_exonends{$l}}=sort { $a <=> $b } @{$transcripts_exonends{$l}};
	@{$transcripts_exonstarts{$l}}=sort { $a <=> $b } @{$transcripts_exonstarts{$l}};
	$start_l=${$transcripts_exonstarts{$l}}[0]-1;

	$end_l=${$transcripts_exonends{$l}}[$#{$transcripts_exonends{$l}}];
	print STDOUT "$transcripts_chr{$l}\t$start_l\t$end_l\t$l\t0\t$transcripts_strand{$l}\n";

}
