#!/usr/bin/env perl

use strict;
use warnings;
$|=1;
use Data::Dumper;

my $sjFile=$ARGV[0];
my $ssFile=$ARGV[1];

#SJ identifier must be in first field of $sjFile.

#print Dumper \%genetypes;
open F, "$sjFile" or die $!;

my %sjLine=();

while (<F>){
	chomp;
	my $line=$_;
	next if ($line=~/^#/);
	my @line=split("\t", $line);
	#die "Wrong number of fields ($#line)\n".join(",", @line)."\n" unless ($#line==4 || $#line==3);
	if($#line>2){
		$line[4]='0' unless (defined $line[4]);
	}
	@{$sjLine{$line[0]}}=@line;
}
#print "SJs\n". Dumper \%sjLine;

open G, "$ssFile" or die $!;
my %minScore=();
$minScore{'Acceptor'}=1000000;
$minScore{'Donor'}=1000000;
#print "minScore B: ".Dumper \%minScore;
my %ssScores=();
while(<G>){
	chomp;
	my $line=$_;
	next if ($line=~/^#/);
	my @line=split("\t", $line);
	my @secondField=split(":", $line[1]);
	my $sj=$secondField[0];
	my $SStype=$line[2];
	my $score=$line[3];
	$ssScores{$sj}{$SStype}=$score;
	if($score<$minScore{$SStype}){
		$minScore{$SStype}=$score;
	}
}

#print "SSscores B\n". Dumper \%ssScores;
#print "minScore A: ".Dumper \%minScore;
print STDERR "Min SS scores found: ";
print STDERR Dumper \%minScore;
print STDERR "\nThese values will be assigned to any SS not found in $ssFile.\n";

my @SStypes=('Acceptor', 'Donor');
#fill in missing values
foreach my $sj (keys %sjLine){
	foreach my $type (@SStypes){
		unless (exists $ssScores{$sj}{$type}){
			$ssScores{$sj}{$type}=$minScore{$type}
		}
	}
}

#print final output with summed scores for each SJ

foreach my $sj (keys %sjLine){
	my $sjScore=$ssScores{$sj}{'Acceptor'}+$ssScores{$sj}{'Donor'};
	print join ("\t", @{$sjLine{$sj}})."\t$sjScore\n";
}
