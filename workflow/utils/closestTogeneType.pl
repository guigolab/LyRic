#!/usr/bin/env perl

use strict;
use warnings;
$|=1;
use Data::Dumper;

my $bedtsvFile=$ARGV[0];
my $genetypesFile=$ARGV[1];
my $clusterIdColumn=4;
$clusterIdColumn=$ARGV[2] if ($ARGV[2]);
print STDERR "clusterIdColumn: $clusterIdColumn\n";
open F, "$genetypesFile" or die $!;

my %genetypes=();
print STDERR "Reading in $genetypesFile...\n";
while(<F>){
	chomp;
	my @line=split "\t";
	my @second=split(",", $line[1]);
	foreach my $sec (@second){
		$sec=~s/.+://g;
		$genetypes{$line[0]}{$sec}=1;
	}
	#print STDERR Dumper \%genetypes;

}
close F;
print STDERR "Done\n";
#print Dumper \%genetypes;
open F, "$bedtsvFile" or die $!;
print STDERR "Reading in $bedtsvFile and writing output...\n";
#my %clusterToRecord=();
#my %lineNumberToRecord=();
my %out=();
while (<F>){
	chomp;
	my $lineNumber=$.;
	my $record=$_;
#	$lineNumberToRecord{$lineNumber}=$record;
	my @line=split "\t";
	my $clusterId=$line[$clusterIdColumn-1]; #$line[3];
	#my $distance=$line[$#line];
	my @readIds=split(/;|,/,$clusterId);
#	if (exists $clusterToRecord{$clusterId}){
#		die "Cluster ID found >1 time in input, cannot continue. Culprit:\n$clusterId\n";
#	}
	my %clusterGeneType=();
	foreach my $read (@readIds){
		if(exists $genetypes{$read}){
			foreach my $genetype (keys %{$genetypes{$read}}){
				$clusterGeneType{$genetype}=1;
				#$clusterToRecord{$clusterId}=$record;
				#$out{$line[3]}{$genetype}=$line[$#line];
				#print "$line[3]\t$genetype\t$line[$#line]\n"
			}

		}
		else{
			print STDERR "$read not found in $genetypesFile, skipped\n";
		}
	}
	my @geneTypes;
	foreach my $key (keys %clusterGeneType){
		push(@geneTypes, $key);
	}
	if($#geneTypes<1){
		print "$record\t$geneTypes[0]\n";
	}
	else{
		print "$record\tmultiBiotype\n";

	}
}
print STDERR "Done\n";

#foreach my $key (keys %clusterToGeneType){
# foreach my $key (sort {$a <=> $b} (keys %lineNumberToRecord)){
# 	my @geneTypes=();
# 	foreach my $key2 (keys %{$clusterToGeneType{$key}}){
# 		push(@geneTypes, $key2);
# 	}
# 	if($#geneTypes<1){
# 		print "$key\t".join(",", @geneTypes)."\t$clusterToRecord{$key}\n";
# 	}
# 	else{
# 		print "$key\tmultiBiotype\t$clusterToRecord{$key}\n";

# 	}
# }
