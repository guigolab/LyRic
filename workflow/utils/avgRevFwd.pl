#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
open FILE1, $ARGV[0] or die "$!";
open FILE2, $ARGV[1] or die "$!";

my %file1=();
while(<FILE1>){
	chomp;
	$_=~s/-nan/0/g;
	my $line=$_;
	my @line=split ("\t", $line);
	my $key=shift(@line);
	#print STDERR "$key\n###\n";
	#print STDERR join("\t", @line);
	@{$file1{$key}}=@line;
}

my %file2=();
while(<FILE2>){
	chomp;
	$_=~s/-nan/0/g;
	my $line=$_;
	my @line=split ("\t", $line);
	my $key=shift(@line);
	#print STDERR "$key\n###\n";
	#print STDERR join("\t", @line);
	@{$file2{$key}}=@line;
}

foreach my $pos (sort { $a <=> $b } keys %file1){
	#print STDERR "$pos\n";
	my @avgs=();

	#for(my $i=0; $i<=$#{$file1{$pos}}; $i++){
	my $avg1=${$file1{$pos}}[0];
	my $avg2=${$file2{$pos}}[0];
	my $size1=${$file1{$pos}}[1];
	my $size2=${$file2{$pos}}[1];
	#artificially set size to pseudocount when it's zero. This will lead to avg=0 instead of pointlessly dividing by zero
	$size1=0.0000001 if $size1==0;
	$size2=0.0000001 if $size2==0;
	my $sum1=$avg1*$size1;
	my $sum2=$avg2*$size2;
	my $avg=($sum1+$sum2)/($size1+$size2);
	push(@avgs, $avg);
	#}
	print "$pos\t".join("\t",@avgs)."\n";
}

