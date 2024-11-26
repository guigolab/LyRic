#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
my $boulderfilename=$ARGV[0];
my $racetype=$ARGV[1];
die "Need to know race type ('5' or '3') as second argument\naborted\n" unless ($racetype eq '5' || $racetype eq '3');
my $strand=undef;
if($racetype eq '5'){
	$strand='-';
}
elsif($racetype eq '3'){
	$strand='+';
}
else{
	die;
}

open PRIMERS, "$boulderfilename" or die $!;
#open GFF, ">$boulderfilename.$racetype.gff" or die $!;
my $targetid='';
my $targetSeq='';

#print STDOUT "#targetId\traceType\tprimerSeq\tprimerRank\tstartOnTarget\tendOnTarget\tprimerTm\tprimerGcContent\n";

my %target2primerSeqs=();
#my %target2primerRanks=();
my %target2primerStarts=();
my %target2primerEnds=();
my %target2primerTm=();
my %target2primerGcPc=();
my %target2primerAmpliconLength=();
my %targetSeq=();
while (<PRIMERS>){
	my $primerrank='';
	my $primerseq='';
	#my $primertype='';
	my $tm='';
	my $size='';
	my $gc_pc='';
	
	if ($_=~/^=$/){
		#end of record
		$targetid=undef;
		$targetSeq=undef;
		next;
	}
	if ($_=~/SEQUENCE_ID=(\S+)\n/){
		$targetid=$1;
	}
	if($_=~/SEQUENCE_TEMPLATE=(\S+)\n/){
		$targetSeq=$1;
		#$targetSeq{$targetid}=$targetSeq;
	}
   	elsif($_=~/PRIMER_INTERNAL_(\d+)_SEQUENCE=(\S+)\n/){
		$primerrank=$1;
		$primerseq=$2;
		if($racetype==5){
			$primerseq=reverse($primerseq);
			$primerseq=~tr/ACGTacgt/TGCAtgca/;
		}
		${$target2primerSeqs{$targetid}}[$primerrank]=$primerseq;
		#print STDOUT ">".$targetid."_primer_$racetype"."race_"."$primerrank\n$primerseq\n";
	}
	
	elsif($_=~/PRIMER_INTERNAL_(\d+)=(\d+),(\d+)\n/){
		my $startontrans=$2;
		my $length=$3;
		$primerrank=$1;
		my $endontrans=$startontrans+$length;
		${$target2primerStarts{$targetid}}[$primerrank]=$startontrans;
		${$target2primerEnds{$targetid}}[$primerrank]=$endontrans;
		if($racetype==5){
			${$target2primerAmpliconLength{$targetid}}[$primerrank]=$endontrans;
		}
		elsif($racetype==3){
			${$target2primerAmpliconLength{$targetid}}[$primerrank]=length($targetSeq)-$startontrans;
		}
		else{
			die;
		}
		#print GFF "$chr\tprimer_$racetype"."prace\tprimer\t$startonchr\t$endonchr\t$primerrank\t.\t.\t".$targetid."_primer_$racetype"."race_"."$primerrank\n";
		
	}
	elsif($_=~/PRIMER_INTERNAL_(\d+)_TM=(\S+)\n/){
		$primerrank=$1;
		${$target2primerTm{$targetid}}[$primerrank]=$2;
	}
	elsif($_=~/PRIMER_INTERNAL_(\d+)_GC_PERCENT=(\S+)\n/){
		$primerrank=$1;
		${$target2primerGcPc{$targetid}}[$primerrank]=$2;
	}
	

}

foreach my $targetid (keys %target2primerSeqs){
	for (my $i=0;$i<=$#{$target2primerSeqs{$targetid}};$i++){
		#print "$targetid\t$racetype\t${$target2primerSeqs{$targetid}}[$i]\t$i\t${$target2primerStarts{$targetid}}[$i]\t${$target2primerEnds{$targetid}}[$i]\t${$target2primerTm{$targetid}}[$i]\t${$target2primerGcPc{$targetid}}[$i]\n";
		my $startGff=${$target2primerStarts{$targetid}}[$i]+1;
		
		print "$targetid\tprimer3\t".$racetype."pRacePrimer\t$startGff\t${$target2primerEnds{$targetid}}[$i]\t.\t$strand\t.\ttranscript_id \"$targetid.$racetype"."p.$i\"; primerSeq \"${$target2primerSeqs{$targetid}}[$i]\"; primerRank \"$i\"; primerTm \"${$target2primerTm{$targetid}}[$i]\"; primerGcContent \"${$target2primerGcPc{$targetid}}[$i]\"; predictedAmpliconLength \"${$target2primerAmpliconLength{$targetid}}[$i]\";\n"; # targetSeq \"$targetSeq{$targetid}\";\n";
	}

}

#print Dumper \%target2primerSeqs;

#print Dumper \%target2primerStarts;
#print Dumper \%target2primerEnds;
