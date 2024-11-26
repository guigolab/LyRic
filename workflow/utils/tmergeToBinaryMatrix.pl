#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Data::Dumper;
$|=1;


my $gtf = $ARGV[0];
my $outPrefix=$ARGV[1];
open F, "$gtf" or die $!;

my %set_to_tids=();
my %tid_to_sets=();
while (<F>){
	chomp;
	$_=~/transcript_id \"(\S+)\";.+contains \"(\S+)\";/;
	my $tid=$1;
	my $ct=$2;
	my @ct=split(',', $ct);
	foreach my $id (@ct){
		$id=~/^=(\S+?)=\S+/;
		my $set=$1;
		$set_to_tids{$set}{$tid}=1;
		$tid_to_sets{$tid}{$set}=1;
	}
}

my %set_sizes=();
foreach my $set (keys %set_to_tids){
	$set_sizes{$set}= keys %{$set_to_tids{$set}};
}


my @sets=sort keys %set_to_tids;
my %intersections=();
for (my $i=0;$i<$#sets;$i++){
	for (my $j=$i+1;$j<=$#sets;$j++){
		$intersections{$sets[$i]}{$sets[$j]}=0;
	}
}
print "transcript_id\t".join ("\t", @sets)."\n";
foreach my $tid (keys %tid_to_sets){
	my @sets_bool=();
	for (my $i=0;$i<=$#sets;$i++){
		my $set=$sets[$i];
		if (exists $tid_to_sets{$tid}{$set}){
			push (@sets_bool, '1');
		}
		else{
			push (@sets_bool, '0')
		}
		for (my $j=$i+1;$j<=$#sets;$j++){
			my $set2=$sets[$j];
			if (exists $tid_to_sets{$tid}{$set} && exists $tid_to_sets{$tid}{$set2}){
				$intersections{$set}{$set2}++;
				$intersections{$set2}{$set}++;
			}
		}
	}
	print "$tid\t".join("\t", @sets_bool)."\n";
	#print Dumper \%intersections;
}
#print Dumper \%intersections;
#print Dumper \%set_to_tids;
open JACCARD, ">$outPrefix.jaccard_ind.tsv";
open OVERLAP, ">$outPrefix.overlap_coeff.tsv";

print JACCARD "\t".join("\t", @sets)."\n";
print OVERLAP "\t".join("\t", @sets)."\n";
for (my $i=0;$i<=$#sets;$i++){
	my @ar=();
	my @jac=();
	my @ov=();
	for (my $j=0;$j<=$#sets;$j++){
		if(exists $intersections{$sets[$i]}{$sets[$j]}){
			push(@ar, $intersections{$sets[$i]}{$sets[$j]});
			#calculate Jaccard index:
			my $jac=$intersections{$sets[$i]}{$sets[$j]} / (($set_sizes{$sets[$i]}+$set_sizes{$sets[$j]})-$intersections{$sets[$i]}{$sets[$j]});
			push(@jac, $jac);
			#calculate overlap coefficient
			my $minSet;
			if($set_sizes{$sets[$i]} <= $set_sizes{$sets[$j]}){
				$minSet = $set_sizes{$sets[$i]};
			}
			else{
				$minSet = $set_sizes{$sets[$j]};
			}
			my $ov_coeff=$intersections{$sets[$i]}{$sets[$j]} / $minSet;
			push(@ov, $ov_coeff);
		}
		elsif ($sets[$i] eq $sets[$j]){
			push(@ar, 1);
			push(@jac, 1);
			push(@ov, 1)

		}
		else{
			push(@ar, 'NA');
			push(@jac, 'NA');
			push(@ov, 'NA')
		}
	}
	print JACCARD "$sets[$i]\t".join("\t",@jac)."\n";
	print OVERLAP "$sets[$i]\t".join("\t",@ov)."\n";
}
