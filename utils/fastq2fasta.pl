#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Data::Dumper;

#open F, "$ARGV[0]" or die;

my $count=0;
my $id;
my $seq;
while (<>){
	chomp;
	my $line=$_;
	$count++;
	if($.%4==1){
		$id=$line;
	}
	elsif($.%4==2){
		$seq=$line;
		print ">$id\n$seq\n"
	}

}
