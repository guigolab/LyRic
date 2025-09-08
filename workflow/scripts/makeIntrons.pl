#!/usr/bin/env perl

use strict;
use warnings;

my $sortedGff;
if(defined($ARGV[0])){
	$sortedGff=$ARGV[0];
}
else{
	die "Need input GTF file name as argument.\n";
}
open GFF, "$sortedGff" or die $!;
my %transcript_exons=();
#my %transcript_chr=();
#my %transcript_strand=();
#my %transcript_gene_id=();
my $previous_start=-1;
my $previous_chr='Caravaggio';

while (<GFF>){
	next if ($_=~/^#/);
	next unless ($_=~/\texon\t/);
	my $line=$_;
	if ($line=~/^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t.*transcript_id "(\S+)?";/){

		if ($3 eq "exon"){
			my $transcript_id=$7;
			my $chr=$1;
			my $start=$4;
			my $stop=$5;
			my $strand=$6;
			my $gene_id=undef;
			unless ($strand eq '+' || $strand eq '-'){
				die "Unrecognized strand value '$strand' at line $.\n";
			}
			die "ERROR: Unsorted GTF input (line $.). Must be sorted by chr, then start, then stop. Can't continue.\n" if ($chr eq $previous_chr && $start < $previous_start);
			if($line=~/^(\S+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t(\S+)\t\S+\t.*gene_id "(\S+)?";/){
				$gene_id=$7;
			}
			#$transcript_strand{$transcript_id}=$strand;
			#$transcript_chr{$transcript_id}=$chr;

			my @exon=($start, $stop);
			push(@{$transcript_exons{$transcript_id}}, \@exon);
			if($#{$transcript_exons{$transcript_id}}>0){
				my $i=$#{$transcript_exons{$transcript_id}}-1;
				my $intronChr=$chr;
				my $intronStrand=$strand;
				my $intronStart=${$transcript_exons{$transcript_id}}[$i][1]+1;
				my $intronStop=${$transcript_exons{$transcript_id}}[$i+1][0]-1;
				my @intron=($intronChr, $intronStart, $intronStop, $intronStrand);
				if(defined $gene_id){
					print "$intron[0]\tmakeIntrons\tintron\t$intron[1]\t$intron[2]\t.\t$intron[3]\t.\tgene_id \"$gene_id\"; transcript_id \"$transcript_id\";\n";
				}
				else{
					print "$intron[0]\tmakeIntrons\tintron\t$intron[1]\t$intron[2]\t.\t$intron[3]\t.\ttranscript_id \"$transcript_id\";\n";
				}
				shift(@{$transcript_exons{$transcript_id}})
			}

			if ($chr ne $previous_chr){
				$previous_chr=$chr;
			}
			$previous_start=$start;
		}
	}
	else{
		die "ERROR: line $.: malformed GTF record. Can't continue.\n";
	}
}
