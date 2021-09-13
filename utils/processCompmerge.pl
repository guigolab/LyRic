#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use FindBin;    # find present script
use lib "$FindBin::Bin";
use gffToHash;
use hashToGff;
use Data::Dumper;

my $compmergeOutGff=$ARGV[0];
my $extraIdprefix='';
$extraIdprefix=$ARGV[1] if (defined ($ARGV[1]));



open GFF, "$compmergeOutGff" or die $!;

my %gffHash=gffToHash($compmergeOutGff, 'transcript_id', 0, 'exon');

#print  Dumper \%gffHash;

foreach my $tx (keys %gffHash){
	foreach my $exon (@{$gffHash{$tx}}){
		#print Dumper $exon;
		#print Dumper @{$exon}[8];
		my $orig_trancript_id=${@{$exon}[8]}{'transcript_id'};
		${@{$exon}[8]}{'transcript_id'}= $extraIdprefix."cm.".$orig_trancript_id.".";
	}

}


my @outGff=hashToGff(\%gffHash);
print join("", @outGff);
