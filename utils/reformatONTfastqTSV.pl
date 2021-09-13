#!/usr/bin/env perl

use strict;

my $fastQtsv=$ARGV[0];
open F, "$fastQtsv" or die $!;

while (<F>){
    my @line=split "\t";
    die "Wrong format" unless $#line==2;
    my @id=split(" ", $line[0]);

    #make up read ID:
    $id[1]=~s/runid=//g;
    $line[0]=join("/", @id[0..1]); 

    #replace Us with Ts in read sequence
    $line[1]=~s/U/T/g;
    print join("\t",@line)

}