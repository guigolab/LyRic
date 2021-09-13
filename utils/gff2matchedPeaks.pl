#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
$|=1;


#arg1=file1
#arg2=file2
#arg3=element id to report (e.g. "transcript_id" value in 9th field)
#arg4=measure to report (e.g. "RPKM" value in 9th field). 

#output is in "matchedPeak" format (tab-separated), to calculate IDR:
# # element_id	value_in_file_1	element_id	value_in_file_2	#peak present in both reps
# # element_id	value_in_file_1 -1	-1				#peak absent in second file
# # -1	-1	element_id	value_in_file_2				#peak absent in first file

open GFF1, "$ARGV[0]" or die $!;
open GFF2, "$ARGV[1]" or die $!;

my $elementAttrid=$ARGV[2];
my $measure=$ARGV[3];
my %elements1=();
my %elements2=();
my %elements1and2=(); #to account for elements absent from either of the two files

while(<GFF1>){
	my@arr=parse_gff($_);
	#print "@arr\n";
	$elements1{$arr[0]}=$arr[1];
	$elements1and2{$arr[0]}=1;
}


while (<GFF2>){
	my @arr=parse_gff($_);
	#print "@arr\n";
	$elements2{$arr[0]}=$arr[1];
	$elements1and2{$arr[0]}=1;
}

#print Dumper \%elements1;

foreach my $el (keys %elements1and2){
	my $el1=-1;
	my $el2=-1;
	my $el1Value=-1;
	my $el2Value=-1;
	if(exists ($elements1{$el})){
		$el1=$el;
		$el1Value=$elements1{$el};
	}
	if(exists ($elements2{$el})){
		$el2=$el;
		$el2Value=$elements2{$el};
	}
	print "$el1\t$el1Value\t$el2\t$el2Value\n";
}



sub parse_gff{
	chomp $_[0];
	$_[0]=~s/\"//g;
	$_[0]=~s/;//g;
	my $elementid=undef;
	my $value=undef;
	if($_[0]=~/\s$elementAttrid (\S+)/){
		$elementid=$1;
		#$elementid=~s/\"//g;
	}
	else{
		die "$elementAttrid attribute not found in line $.. Cannot continue\n";
	}
	if ($_[0]=~/\s$measure (\S+)/){
		$value=$1;
#		unless($value =~ /^-?\d+\.?\d*$/){
#			warn "Value '$value' is not a number at line $. This will create problems during IDR calculations downstream.";
#		}
	}
	else{
		die "$measure attribute not found in line $.. Cannot continue\n";
	}
	return ($elementid, $value)
}
