#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;    # find present script
use lib "$FindBin::Bin";
use gffToHash;
use hashToGff;
my %gffHash=gffToHash($ARGV[0], 'transcript_id', 0,);

my @outGff=hashToGff(\%gffHash);
print join("", @outGff);

