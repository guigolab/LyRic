#!/usr/bin/env perl

@keys=@ARGV;
#print STDERR "@keys\n";

while(<STDIN>){
	@line=split "\t";
	@out=();
	foreach $key (@keys){

		if($line[8]=~/(\s$key|^$key) "(\S+)";*/ || $line[8]=~/(\s$key|^$key) (\S+);*/){
			push(@out,$2);
		}
		else{
			warn "key $key not found:\n$_"
		}
	}
	if($#out>=0){
		print join("\t", @out)."\n";
	}
}
