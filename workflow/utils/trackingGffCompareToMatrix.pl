#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
$|=1;


my $cmd=$ARGV[0];
my $tracking=$ARGV[1];

open F, "$cmd" or die $!;
my @inputGtfs=();
while (<F>){
	chomp;
	next if($_=~/Command line was:/);
	if($_=~/#gffcompare (.+)/){
		my @args=split(/\s+/, $1);
		for (my $i=0; $i<=$#args; $i++){
			if($args[$i]=~/^-/){ #arg is an option, skip
				$i++
			}
			else{
				push(@inputGtfs, $args[$i]);
			}
		}
		last;
	}
}

close F;
print STDERR scalar @inputGtfs." input files according to $cmd\n";

open F, "$tracking" or die $!;

#my @combTranscripts=();
print STDOUT "transcript_id\t".join("\t",@inputGtfs)."\n";
while (<F>){
	chomp;
	my @combTranscripts=split "\t";
	my $id=$combTranscripts[0];
	for (my $i=0;$i<=3;$i++){
		shift @combTranscripts;
	}
	if(scalar @inputGtfs != scalar @combTranscripts){
		die "Number of input files do not match, cannot continue.\n";
	}
	#print scalar @combTranscripts." elements\n";
	my @boolCombTranscripts=();
	foreach my $file (@combTranscripts){
		if($file eq '-'){
			push(@boolCombTranscripts, '0')
		}
		else{
			push(@boolCombTranscripts, '1')
		}
	}
	print STDOUT $id."\t".join("\t",@boolCombTranscripts)."\n";
}
