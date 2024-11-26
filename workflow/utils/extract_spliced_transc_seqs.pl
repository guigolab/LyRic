#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Bio::DB::Fasta;

# 2 args:
# arg1 -> GFF file containing transcripts info (i.e. exon records, grouped by transcript_id)
# arg2 -> multi-fasta file containing all chromosomes present in the GFF

open I, "$ARGV[0]" or die $!;

my $genomeFa=$ARGV[1];
my $chrdb = Bio::DB::Fasta->new("$genomeFa");
my %transcripts_exonstarts;
my %transcripts_exonends;
my %transcripts_strand;
my %transcripts_chr;

print STDERR "Reading GFF in ($ARGV[0])...\n";
while (<I>){
	if ($_=~/^(\S+)\t(\S+)\texon\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*transcript_id "(\S+)?";.*\n/ || $_=~/^(\S+)\t(\S+)\texon\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t(\S+).*\n/){
		my $chr=$1;
		my $source=$2;
		my $start=$3;
		my $end=$4;
		my $strand=$5;
		my $transc=$7;
		push(@{$transcripts_exonstarts{$transc}},$start);
		push(@{$transcripts_exonends{$transc}},$end);
		$transcripts_strand{$transc}=$strand;
		$transcripts_chr{$transc}=$chr;

	}
}
my $countTranscripts=0;
foreach my $tr (keys %transcripts_chr){
	$countTranscripts++;
}
print STDERR "Done ($countTranscripts transcripts found in GFF file).\nNow extracting sequences...\n";

my $countTranscripts2=0;
foreach my $transcript (keys %transcripts_exonstarts){
	$countTranscripts2++;
	my $seq='';
	my $chr=$transcripts_chr{$transcript};
	unless($#{$transcripts_exonstarts{$transcript}}==$#{$transcripts_exonends{$transcript}}){
		print STDERR "Nr of exon starts differs from nr of exon ends for transcript $transcript!!!\n";
	}
	if ($transcripts_strand{$transcript} eq '+' || $transcripts_strand{$transcript} eq '.' ){
		@{$transcripts_exonends{$transcript}}=sort { $a <=> $b } @{$transcripts_exonends{$transcript}};
		@{$transcripts_exonstarts{$transcript}}=sort { $a <=> $b } @{$transcripts_exonstarts{$transcript}};

	}
	elsif ($transcripts_strand{$transcript} eq '-'){
		@{$transcripts_exonends{$transcript}}=sort { $b <=> $a } @{$transcripts_exonends{$transcript}};
		@{$transcripts_exonstarts{$transcript}}=sort { $b <=> $a } @{$transcripts_exonstarts{$transcript}};

	}
	for (my $i=0;$i<=$#{$transcripts_exonstarts{$transcript}};$i++){
		my $start=${$transcripts_exonstarts{$transcript}}[$i];
		my $end=${$transcripts_exonends{$transcript}}[$i];
		if($transcripts_strand{$transcript} eq '+'){
			my $tmp=$chrdb->seq($chr, $start => $end);
			if(defined ($tmp)){
				$seq.=$tmp;
			}
			else{
				warn "$transcript $chr $start $end: retrieved undefined sequence\n";
			}
		}
		elsif($transcripts_strand{$transcript} eq '.' ){
			warn "Transcript '$transcript' has strand '.'. Strand converted to '+' for seq extraction.\n";
			my $tmp=$chrdb->seq($chr, $start => $end);
			if(defined ($tmp)){
				$seq.=$tmp;
			}
			else{
				warn "$transcript $chr $start $end: retrieved undefined sequence\n";
			}
		}
		elsif($transcripts_strand{$transcript} eq '-'){
			my $tmp=$chrdb->seq($chr, $start => $end);
			if(defined ($tmp)){
				my $rcseq.=$tmp;
				my $rseq=reverse($rcseq);
				$rseq=~ tr/ACGTacgt/TGCAtgca/;
				$seq.=$rseq;
			}
			else{
				warn "$transcript $chr $start $end: retrieved undefined sequence\n";
			}
		}
	}
	$seq=~s/\s//g;
	$seq=uc($seq);
	print ">$transcript\n";
	for(my $k=0;$k<=length($seq);$k=$k+60){
		print substr($seq,$k,60)."\n";
	}
	if ($countTranscripts2%1000 == 0){
		print STDERR "\tProcessed $countTranscripts2 transcripts\n";
	}
}

print STDERR "Done.\n";