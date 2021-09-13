#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

open TSV, "$ARGV[1]" or die $!;

my %transcriptToRightStrand=();

while(<TSV>){
	chomp;
	next if ($_=~/^#/);
	my @line=split "\t";
	if (exists $transcriptToRightStrand{$line[0]}){
		die "$line[0] present at least twice. Can't continue.\n";
	}
	$transcriptToRightStrand{$line[0]}=$line[1]
}
close TSV;
open GFF, "$ARGV[0]" or die $!;
#my $notFound=$ARGV[0].".notfound.gff";
#open NOTFOUND, ">$notFound" or die $!;
while(<GFF>){
	next if $_=~/^#/;
	my @line=split "\t";
	$line[8]=~/transcript_id \"(\S+)\"/;
	my $transcript_id=$1;
	if(exists $transcriptToRightStrand{$transcript_id}){
		$line[6]=$transcriptToRightStrand{$transcript_id};
		print join ("\t", @line);

	}
#	else{
#		print NOTFOUND join ("\t", @line);

#	}
}
close GFF;
#close NOTFOUND;
