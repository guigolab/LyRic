#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
#use re 'debug';
$/=undef;

`mkdir -p frames/`;
`rm frames/*`;

my $line=<>;
my %frames=();
#if($_=~/\\begin\{frame\}(.+)\\end\{frame\}/){
if (my @frames = $line =~ /\\begin\{frame\}.*?\\end\{frame\}/sg) {
	foreach my $frame (@frames){
		my $title='dummy';
		if($frame =~ /\\frametitle\{(.*?)\}\n/){
			$title=$1;
			#print "###\n$title\n";
			$title=~s/[^A-Za-z0-9]//g;
		}
		#print "$title\n";
		$title=newKeyIfExists($title);
		$frames{$title}=$frame;

	}
	#print join("\n###\n", @frames)
	#print "\n\n####### BEGIN ########\n\n$1\n\n####### END #########\n\n"
}

foreach my $key (keys %frames){
	print STDERR "Frame: $key\n";
	my $outfile="frames/$key.tex";
	open OUT, ">$outfile" or die $!;
	print OUT $frames{$key};
	close OUT;
}

sub newKeyIfExists{
	my $key=$_[0];
	#print "key $key\n";
	if(exists $frames{$key}){
		$key.=".1";
		newKeyIfExists($key);
	}
	else{
		return $key;
	}
}