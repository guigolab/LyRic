#!/usr/bin/env perl
use strict;
use warnings;
use lib "/users/rg/jlagarde/julien_utils/";
use gffToHash;
use hashToGff;



my $gffIn=$ARGV[0];


print STDERR "Parsing input GFF $gffIn\n";
my %gffHash=gffToHash($gffIn, 'transcript_id', 1, 'exon');
print STDERR "Done\n";
#print STDERR Dumper \%gffHash;

#my %outGff=dclone(%gffHash); #we'll only edit outGff

#print STDERR Dumper \%gffHash;
my @outGff=hashToGff(\%gffHash);
print join("", @outGff);
