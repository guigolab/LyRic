#!/usr/bin/env perl

use strict;
use warnings;

my $bed=$ARGV[0];
open F, "$bed" or die $!;

while (<F>){
	my @line=split "\t";
	my @blockSizes=split(",",$line[10]);
	my $splicedLength=0;
	foreach my $block (@blockSizes){
		$splicedLength+=$block

	}
	print STDOUT "$line[3]\t$splicedLength\n";
}
