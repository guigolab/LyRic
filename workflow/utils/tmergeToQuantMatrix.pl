#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Data::Dumper;
$|=1;


my $gtf = $ARGV[0];
my $quantTsv = $ARGV[1];

open Q, "$quantTsv" or die $!;
my %tidToExpValue=();

while (<Q>){
    chomp;
    my @line=split("\t");
    $tidToExpValue{$line[0]}=$line[1];
}
close Q;

open F, "$gtf" or die $!;

my %tid_to_sets=();
my %sets=();
while (<F>){
	chomp;
	$_=~/transcript_id \"(\S+)\";.+contains \"(\S+)\";/;
	my $tid=$1;
	my $ct=$2;
	my @ct=split(',', $ct);
	foreach my $id (@ct){
		$id=~/^=(\S+?)=\S+/;
		my $set=$1;
        $sets{$set}=undef;
        die "Could not find expression value for $id in $quantTsv\n" unless (exists $tidToExpValue{$id});
        my $expValue=$tidToExpValue{$id};
		$tid_to_sets{$tid}{$set}=$expValue;
	}
}
%tidToExpValue=();
my @sets=sort keys %sets;


print "transcript_id\t".join ("\t", @sets)."\n";
foreach my $tid (keys %tid_to_sets){
	my @sets_expValue=();
	for (my $i=0;$i<=$#sets;$i++){
		my $set=$sets[$i];
		if (exists $tid_to_sets{$tid}{$set}){
			push (@sets_expValue, $tid_to_sets{$tid}{$set});
		}
		else{
			push (@sets_expValue, '0')
		}
	}
	print "$tid\t".join("\t", @sets_expValue)."\n";
	#print Dumper \%intersections;
}
