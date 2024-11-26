#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my $binSize=1;
my $numberOfBins=100;

open TMERGE, "$ARGV[0]" or die $!;

my %transcript_to_threepDists=();
my %transcript_to_fivepDists=();
while (<TMERGE>){
    #print;
    chomp;
    $_=~/transcript_id "(\S+)";/;
    my $trId=$1;
    unless(exists $transcript_to_threepDists{$trId}){
        $_=~/meta_3p_dists_to_5p "(\S+)";/;
        @{$transcript_to_threepDists{$trId}}=split(",", $1);
    }
    unless(exists $transcript_to_fivepDists{$trId}){
        $_=~/meta_5p_dists_to_5p "(\S+)";/;
        @{$transcript_to_fivepDists{$trId}}=split(",", $1);
    }
}

close TMERGE;
#print Dumper \%transcript_to_fivepDists;

my %metaGeneBins=();
for (my $i=0; $i<$numberOfBins; $i=$i+$binSize){
    $metaGeneBins{$i}=0;
}
print "bin\tcount\n";
foreach my $trId (keys %transcript_to_fivepDists){
    for (my $i=0; $i<= $#{$transcript_to_fivepDists{$trId}}; $i++){
        my $fivepDist=$numberOfBins*${$transcript_to_fivepDists{$trId}}[$i];
        my $threepDist=$numberOfBins*${$transcript_to_threepDists{$trId}}[$i];
        #print STDERR "$trId\t$fivepDist\t$threepDist\n";
        foreach my $bin (sort { $a <=> $b } keys %metaGeneBins){
            if($fivepDist < $bin + $binSize && $threepDist >= $bin){
                $metaGeneBins{$bin}++;
            }
        }
    }

}

foreach my $bin (sort { $a <=> $b } keys %metaGeneBins){
    print "$bin\t$metaGeneBins{$bin}\n";
}