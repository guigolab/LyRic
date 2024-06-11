#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use lib "/users/rg/jlagarde/julien_utils/";
use reverseComplement;

my $motif=uc($ARGV[0]);
my $revCompMotif=reverseComplement($motif);
open TBL, "$ARGV[1]" or die $!;
#takes tbl as input

while (<TBL>){
	chomp;
	my @line=split " ";
	my $chr=$line[0];
	print STDERR "Chr: $line[0]...\n";
	my $seq=uc($line[1]);
	while($seq =~/$motif/g) {
		print "$chr\t".($-[0])."\t".($+[0])."\t".$motif.".".$chr."_".($-[0])."_".($+[0])."_+\t0\t+\n";
    }
    while($seq =~/$revCompMotif/g) {
		print "$chr\t".($-[0])."\t".($+[0])."\t".$motif.".".$chr."_".($-[0])."_".($+[0])."_-\t0\t-\n";
    }
}