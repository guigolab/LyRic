#!/usr/bin/env perl

### When replicates are merged by tmerge, e.g.:

# gtf_sampleA_rep1 ----tmerge----> gtf_sampleA_rep1_merged

# gtf_sampleA_rep2 ----tmerge----> gtf_sampleA_rep2_merged
# (THIS IS STEP 1)

### One may want to further merge these output:

# first concatenate replicate files:
# cat gtf_sampleA_rep2_merged gtf_sampleA_rep1_merged |sortgff > gtf_sampleA_rep1_rep2_merged

# then merge:
# tmerge gtf_sampleA_rep1_rep2_merged > gtf_sampleA_merged
# (THIS IS STEP 2)

# The 'contains' GTF attribute of TMs represented in gtf_sampleA_merged includes a list of transcript_id's from gtf_sampleA_rep1_rep2_merged, which is not always desirable. 

# This script recovers transcript_ids listed in the 'contains' attribute of each TM of gtf_sampleA_rep1_rep2_merged and associates them to transcript_id's of gtf_sampleA_merged in a tab-separated file.

# Note that step one can involve a theoretically infinite number of replicates to be merged in step 2.

use strict;
use warnings FATAL => 'all';
use Data::Dumper;
$|=1;


my $gtf1 = $ARGV[0];
my $gtf2 = $ARGV[1];
open F1, "$gtf1" or die $!;

my %gtf1TmToContains=();

while (<F1>){
	chomp;
	$_=~/transcript_id \"(\S+)\";.+contains \"(\S+)\";/;
	my $tid=$1;
	my $ct=$2;
	my @ct=split(',', $ct);
	$gtf1TmToContains{$tid}=\@ct;
}

close F1;

#print Dumper \%gtf1TmToContains;


open F2, "$gtf2" or die $!;

my %gtf2TmToContains=();

while (<F2>){
	chomp;
	$_=~/transcript_id \"(\S+)\";.+contains \"(\S+)\";/;
	my $tid=$1;
	my $ct=$2;
	my @ct=split(',', $ct);
	$gtf2TmToContains{$tid}=\@ct;
}

close F2;
foreach my $tid3 (keys %gtf2TmToContains){
	foreach my $tid2 (@{$gtf2TmToContains{$tid3}}){
		unless (exists $gtf1TmToContains{$tid2}){
			die "Couldn't find transcript_id $tid2 in $gtf1.\n";
		}
		foreach my $tid1 (@{$gtf1TmToContains{$tid2}}){
			print "$tid1\t$tid3\n";
		}
	}
}

