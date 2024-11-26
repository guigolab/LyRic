#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;



my $color_mappings=$ARGV[0];
my $bedFile=$ARGV[1];
my $comptr_out=$ARGV[2];
my $trackName=basename($bedFile);

open F, "$color_mappings" or die $!;
my %colors=();
while (<F>){
	chomp;
	my @line=split "\t";
	$colors{$line[0]}=$line[1];
}
close F;



open COMPTR, "$comptr_out" or die $!;
my %tr_to_comptr=();
while(<COMPTR>){
	chomp;
	my @line=split "\t";
	$line[0]=~s/"//g;
	$line[0]=~s/;//g;
	$tr_to_comptr{$line[0]}=$line[1];
}

open BED, "$bedFile" or die $!;

print STDOUT "track name=$trackName description=\"$trackName\" itemRgb=On\n";

while(<BED>){
	chomp;
	my @line=split "\t";
	my $trId=$line[3];
	my $comptrCategory=$tr_to_comptr{$trId};
	my $rgb=$colors{$comptrCategory};
	$line[8]=$rgb;
	print STDOUT join("\t", @line)."\n";
}
