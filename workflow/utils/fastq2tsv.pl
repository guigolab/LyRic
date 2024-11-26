#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;
use Data::Dumper;

#open F, "$ARGV[0]" or die;

my $count=0;
my $id;
my $seq;
my $qual;
while (<>){
	chomp;
	my $line=$_;
	$count++;
	if($.%4==1){ #read id
		$id=$line;
		$id=~s/^@//;
	}
	elsif($.%4==2){ # read seq
		$seq=$line;
	}
	elsif($.%4==0){ # qualities
		$qual=$line;
		print "$id\t$seq\t$qual\n"
	}

}
