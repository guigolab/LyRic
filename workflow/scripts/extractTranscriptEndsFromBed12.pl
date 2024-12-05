#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
die "Specify 5 or 3 (ie 3' or 5' end to extract).\n" unless ($ARGV[0] && ($ARGV[0] == 3 || $ARGV[0] == 5));

my $endToExtract=$ARGV[0];


while (<STDIN>){
	chomp;
	my @line=split "\t";
	if($endToExtract == 5){
		if ($line[5] eq '+'){
			$line[2]=$line[1]+1;
		}
		elsif($line[5] eq '-'){
			$line[1]=$line[2]-1;
			$line[2]=$line[2];
		}
		else{
			die "Line $. of input: strand (field 6) must be '+' or '-'. Cannot continue.\n";
		}
	}
	else{
		if($line[5] eq '-'){
			
			$line[2]=$line[1]+1;
		}
		elsif($line[5] eq '+'){
			$line[1]=$line[2]-1;
			$line[2]=$line[2];
		}
		else{
			die "Line $. of input: strand (field 6) must be '+' or '-'. Cannot continue.\n";
		}
	}
	splice (@line, 6);
	print join ("\t", @line)."\n";
}
