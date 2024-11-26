#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
$|=1;

#input must be BLAST tabular output  with option: '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen score gaps" '

my $minBarcodeHspCoverage=0.7; # WARNING: if barcode is in the middle of the adapters, be sure to require >> 50% of coverage, otherwise many duplicates will be found at a given position!
my $evalue_cutoff=1;

my $count_discarded_hits=0;

open F, "$ARGV[0]" or die $!;

my %hits=();
my %read_length=();

while (<F>){
	next if $_=~/^#/;
	chomp;
	my @line=split "\t";
	my $query_id=$line[0]; # barcode id
	my $subject_id = $line[1]; # read id
	my $aln_length = $line[3]; # alignment length
	my $query_length=$line[12];
	my $evalue=$line[10];
	# filter insignificant hits based on % of barcode sequence aligned and e-value
	if($aln_length/$query_length<=$minBarcodeHspCoverage || $evalue>=$evalue_cutoff){
		$count_discarded_hits++;
		next;
	}
	my $aln_score = $line[11];
	my $subject_start = $line[8];
	my $subject_end = $line[9];
	if($subject_end<$subject_start) { # if match is on minus strand of subject (aka read, swap start and end
		my $tmp=$subject_start;
		$subject_start=$subject_end;
		$subject_end=$tmp;
	}
	my $matchMidCoordinate= ($subject_start + $subject_end) /2;
	my $subject_length = $line[13];
	$read_length{$subject_id}=$subject_length;
	my $normalized_location=sprintf("%.2f",$matchMidCoordinate / $subject_length);
	#my @summaryAln = ($aln_score, $subject_start, $subject_end);
	#print join("\t", @line)."\t$normalized_location\n";
	$hits{$subject_id}{$normalized_location}{$query_id} = $aln_score;

}


#print Dumper \%hits;

foreach my $read (keys %hits){
	foreach my $location (keys $hits{$read}){
		my $previous_score= -10;
		my $output_line;
		foreach my $barcode ( sort { $hits{$read}{$location}{$b} <=> $hits{$read}{$location}{$a} } keys %{$hits{$read}{$location}}){
			my $score = $hits{$read}{$location}{$barcode};
			if($score < $previous_score){

				last;
			}
			elsif($score == $previous_score){
				print STDERR "Too many ties at location $location on read $read. Skipped.\n";
				$output_line=undef;
				last;
			}
			else{
				$output_line= "$read\t$location\t$barcode\t$score\n";
				$previous_score=$score;
			}
		}
		if(defined ($output_line)){
			print $output_line;
		}
	}

}

if($count_discarded_hits>0){
	print STDERR "
##################################################################
##################################################################
########## $count_discarded_hits hits were discarded as a result of
########## the e-value and minimum alignment coverage cutoff.
########## These cutoffs are:
##########  - e-value must be < $evalue_cutoff
##########  - aligned portion of barcode sequence must be >= $minBarcodeHspCoverage
##################################################################
##################################################################\n";
}


print STDERR "Done.\n";