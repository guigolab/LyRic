#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Bio::DB::Fasta;



# 2 args:
# arg1 -> BED6 file containing sequence interval
# arg2 -> multi-fasta file containing all chromosomes present in the BED

#extract genomic sequence of each record of exon BED file.
#prints FASTA to stdout


open I, "$ARGV[0]" or die $!;

die "first argument: bed6 file, second argument: chr mmultifasta\naborted\n" unless $ARGV[1];
my $genomeFa=$ARGV[1];
my $chrdb = Bio::DB::Fasta->new("$genomeFa");


#open NUC, ">$nuc_out" or die $!;

while (<I>){
	chomp;
	my @line=split "\t";
	my $chr=$line[0];
	my $start=$line[1]+1;
	my $end=$line[2];
	my $name=$line[3];
	my $score=$line[4];
	my $strand=$line[5];
	my $seq='';
	if($strand eq '+'){
		$seq=$chrdb->seq($chr, $start => $end);
	}
	elsif($strand eq '-'){
		my $rcseq.=$chrdb->seq($chr, $start => $end);
		my $rseq=reverse($rcseq);
		$rseq=~ tr/ACGTacgt/TGCAtgca/;
		$seq.=$rseq;
	}
	else{
		die "Strand must be '+' or '-'. Cannot continue.\n"
	}
	$seq=~s/\s//g;
	$seq=uc($seq);
	print ">$name\n";
	for(my $k=0;$k<=length($seq);$k=$k+60){
		print substr($seq,$k,60)."\n";
	}
}
