#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;
use Data::Dumper;

#open F, "$ARGV[0]" or die;

my $id;
my $seq;
my $qual;
while (<>){
	chomp;
	my @line=split "\t";
	my $id=$line[0];
	unless($id=~/^@/){
		my $tmp="@".$id;
		$id=$tmp;
	}
	my $seq=$line[1];
	my $qual='';
	if ($line[2]){
		$qual=$line[2];
	}
	else{ #if quality scores are absent, make up string with maximum score (40, i.e. 'I') all along
		my $numOfChar = length($seq);
		for (my $i=0;$i<$numOfChar;$i++){
			$qual.='I'
		}
	}
	print "$id
$seq
+
$qual
";

}
