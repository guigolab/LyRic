#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use diagnostics;
$|=1;

#removes terminal introns from an intron list formatted like:
# BG393879        chr20   33398490        33398650        160
# (only the 3 first fields are mandatory)

open LIST, "$ARGV[0]" or die $!;

my %hit_intron_starts;
my %intron_line;
while (<LIST>){
	my @line=split "\t";
	my $line=$_;
	push(@{$hit_intron_starts{$line[0]}}, $line[2]);
	$intron_line{$line[0]}{$line[2]}=$line;
}
close LIST;
foreach my $hit (keys %hit_intron_starts){
	#print STDERR "$hit\n";
	@{$hit_intron_starts{$hit}}= sort { $a <=> $b } (@{$hit_intron_starts{$hit}});
	for(my $i=1; $i<$#{$hit_intron_starts{$hit}};$i++){ #we want to skip first and last intron of each hit
		print STDOUT $intron_line{$hit}{${$hit_intron_starts{$hit}}[$i]};
	}
}
