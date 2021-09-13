#!/usr/bin/perl -w

use warnings;
use strict;
use Data::Dumper;

#input sorted GTF!!!

my %transcripts_exon_starts=();
my %transcripts_exon_ends=();
my %transcript_chr=();
my %transcript_cat=();
my %transcript_str=();
my %transcript_length=();

while(<>){
	chomp;
	my @line=split "\t";
	if($line[2] eq 'exon'){
		$line[8]=~/transcript_id "(\S+)";/;
		my $transcript_id=$1;
		unless(exists $transcript_length{$transcript_id}){
			$transcript_length{$transcript_id}=0;
		}
		push(@{$transcripts_exon_starts{$transcript_id}}, $line[3]);
		push(@{$transcripts_exon_ends{$transcript_id}}, $line[4]);
		$transcript_chr{$transcript_id}=$line[0];
		$transcript_cat{$transcript_id}=$line[1];
		$transcript_str{$transcript_id}=$line[6];
		$transcript_length{$transcript_id}+=($line[4]-$line[3])+1;
	}
	
}
print "#transcript_id\tIntronCoordOnTranscript\tIntronLengthOnGenome\n";
foreach my $tr (keys %transcript_chr){
	#skip monoexonic
	if ($#{$transcripts_exon_starts{$tr}}<1){
		print STDERR "Skipped $tr, as it doesn't contain any intron\n";
		next;
	}
	
	my $current_coord=0;
	if($transcript_str{$tr} eq '+' || $transcript_str{$tr} eq '.'){
		for(my $i=0;$i<$#{$transcripts_exon_starts{$tr}};$i++){
			my $intron_length=(${$transcripts_exon_starts{$tr}}[$i+1]-${$transcripts_exon_ends{$tr}}[$i])+1;
			#my $exonlengthminus1=${$transcripts_exon_ends{$tr}}[$i]-${$transcripts_exon_starts{$tr}}[$i];
			my $exonlength=(${$transcripts_exon_ends{$tr}}[$i]-${$transcripts_exon_starts{$tr}}[$i])+1;
			$current_coord+=$exonlength;
			print "$tr\t$current_coord\t$intron_length\n";
			
		}
	}
	else{
		for(my $i=$#{$transcripts_exon_starts{$tr}};$i>0;$i--){
			my $intron_length=(${$transcripts_exon_starts{$tr}}[$i]-${$transcripts_exon_ends{$tr}}[$i-1])+1;
			#my $exonlengthminus1=${$transcripts_exon_ends{$tr}}[$i]-${$transcripts_exon_starts{$tr}}[$i];
			my $exonlength=(${$transcripts_exon_ends{$tr}}[$i]-${$transcripts_exon_starts{$tr}}[$i])+1;
			$current_coord+=$exonlength;
			print "$tr\t$current_coord\t$intron_length\n";
		}
	}
	

}
