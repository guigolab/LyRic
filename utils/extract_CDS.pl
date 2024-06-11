#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
#use Bio::DB::GFF;
#use Bio::SeqIO;
use Data::Dumper;

## given a gtf file this script outputs AA sequences as well as spliced nucl. sequences of all CDS.
## skips 'pseudogene' records
## 


#print STDOUT "Translate using provided frame info for CDS start when start_codon is missing? [y/n]\n";

my $frameinfo=$ARGV[1];
chomp $frameinfo;
unless($frameinfo eq 'y' || $frameinfo eq 'n'){
	die;
}
my $seqdir=$ARGV[2];



open I, "$ARGV[0]" or die $!;

my %transcripts_cdsstarts;
my %transcripts_cdsends;
my %transcripts_strand;
my %transcripts_chr;
my %transc_stop_codon;
my %transc_start_codon;
my $nuc_out=$ARGV[0]."_CDS.fa";
my $pep_out=$ARGV[0]."_CDS.pep";
my $readme_out=$ARGV[0]."_CDS.readme";
my %frame;
open NUC, ">$nuc_out" or die $!;
open PEP, ">$pep_out" or die $!;
open README, ">$readme_out" or die $!;
print README 
"$nuc_out contains spliced nucleotide sequences of all CDS records contained in file $ARGV[0],
 in fasta format (hg17/NCBI 35 build).
 When a start codon/stop codon was not explicitly annotated 
(i.e. absence of a 'start_codon'/'stop_codon' record in source gff file 
 for a particular transcript), the '[no_start]'/'[no_stop]' flag 
is added to the fasta header of this transcript.\n\n";
print README
"$pep_out contains AA sequences of all CDS records contained in file $ARGV[0], 
in fasta format (hg17/NCBI 35 build).
 Stop codon is not included.
 When a start codon codon was not explicitly annotated 
 (i.e. absence of a 'start_codon' record in source gff file 
 for a particular transcript), ";
if ($frameinfo eq 'n'){
 print README "the spliced nucleotide sequence of the CDS is 
 translated in 3 frames, and only the longest AA sequence is kept. Additionally, ";
}
print README " the '[no_start]' flag is added to the fasta header of this transcript.
 When a stop codon was not explicitly annotated (i.e. absence of a 'stop_codon' record in source gff file 
 for a particular transcript), the '[no_stop]' flag is added to the fasta header of this transcript.\n";

while (<I>){
	if ($_=~/^(\S+)\t(\S+)\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*transcript_id "(\S+)?";.*\n/ || $_=~/^(\S+)\t(\S+)\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t(\S+).*\n/){
		my $chr=$1;
		my $source=$2;
		unless($source=~/pseudogene/){
			my $start=$3;
			my $end=$4;
			my $strand=$5;
			my $frame=$6;
			my $transc=$7;
			#$vega2gbrowse{$transc}=$source."_".$transc;
			#push(@transcs, $transc);
			push(@{$transcripts_cdsstarts{$transc}},$start);
			push(@{$transcripts_cdsends{$transc}},$end);
			$transcripts_strand{$transc}=$strand;
			$transcripts_chr{$transc}=$chr;
			${$frame{$transc}}{$start}=$frame;
		}
		
	}
	elsif($_=~/^(\S+)\t(\S+)\tstop_codon\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*transcript_id "(\S+)?";.*\n/){
		push(@{$transc_stop_codon{$7}}, ($3, $4));
	}
	elsif($_=~/^(\S+)\t(\S+)\tstart_codon\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*transcript_id "(\S+)?";.*\n/){
		push(@{$transc_start_codon{$7}}, ($3, $4));
	}

}


foreach my $transcript (keys %transcripts_cdsstarts){
	my $fasta_remark='';
	print STDOUT "$transcript\n";
	unless(exists $transc_start_codon{$transcript}){
		$fasta_remark.=' [no_start]';
		#print STDERR "\tWARNING: transcript $transcript 's  start codon is not explicitly annotated, translating anyway\n";
		
	}

	unless(exists $transc_stop_codon{$transcript}){
		#print STDERR "\tWARNING: transcript $transcript 's  stop codon is not explicitly annotated, translating anyway\n";
		$fasta_remark.=' [no_stop]';
	}
	my $seq='';
	my $pepseq='';
	my $chr_file=$transcripts_chr{$transcript}.".fa";
	my $first_frame;
	unless($#{$transcripts_cdsstarts{$transcript}}==$#{$transcripts_cdsends{$transcript}}){
		print STDERR "Nr of CDS starts differs from nr of CDS ends for transcript $transcript!!!\n";
	}
	if ($transcripts_strand{$transcript} eq '+'){
		@{$transcripts_cdsends{$transcript}}=sort { $a <=> $b } @{$transcripts_cdsends{$transcript}};
		@{$transcripts_cdsstarts{$transcript}}=sort { $a <=> $b } @{$transcripts_cdsstarts{$transcript}};
		
	}
	elsif ($transcripts_strand{$transcript} eq '-'){
		@{$transcripts_cdsends{$transcript}}=sort { $b <=> $a } @{$transcripts_cdsends{$transcript}};
		@{$transcripts_cdsstarts{$transcript}}=sort { $b <=> $a } @{$transcripts_cdsstarts{$transcript}};
		
	}
	
	else{
		print STDERR "Couldn't find strand info for transcript $transcript\n";
	}
	$first_frame=${$frame{$transcript}}{${$transcripts_cdsstarts{$transcript}}[0]}+1;
	for (my $i=0;$i<=$#{$transcripts_cdsstarts{$transcript}};$i++){
		
				
		my $start=${$transcripts_cdsstarts{$transcript}}[$i];
		my $end=${$transcripts_cdsends{$transcript}}[$i];
		
#		if ($i==$#{$transcripts_cdsstarts{$transcript}}){ #to include stop codon in output seq
#			if($transcripts_strand{$transcript} eq '+'){
#				if(exists @{$transc_stop_codon{$transcript}}){
					
#				$end=$end+3;
#			}
#			elsif($transcripts_strand{$transcript} eq '-'){
#				$start=$start-3;
#			}
#		}


		if($transcripts_strand{$transcript} eq '+'){
			$seq.= `chr_subseq $seqdir/$chr_file $start $end |extractseq -auto -filter -osformat2 plain `;
		}
		elsif($transcripts_strand{$transcript} eq '-'){
			#print STDERR 
			$seq.=`chr_subseq $seqdir/$chr_file $start $end | revseq -auto -filter -osformat2 plain`;
		}
	}
	
	$seq=~s/\s//g;
	print NUC ">$transcript$fasta_remark\n";
	
	for(my $k=0;$k<=length($seq);$k=$k+60){
		print NUC substr($seq,$k,60)."\n";
	}
	my $pep;
my $frame_used='';
	if ($frameinfo eq 'n' && !exists $transc_start_codon{$transcript}){
		my $pep1=`transeq -frame 1 -auto -stdout -osformat2 plain $nuc_out:$transcript`;
		my $pep2=`transeq -frame 2 -auto -stdout -osformat2 plain $nuc_out:$transcript`;
		my $pep3=`transeq -frame 3 -auto -stdout -osformat2 plain $nuc_out:$transcript`;
		my@pep1 = split(//, $pep1);
		my@pep2 = split(//, $pep2);
		my@pep3 = split(//, $pep3);
		my $stop_idx1=$#pep1;
		my $stop_idx2=$#pep2;
		my $stop_idx3=$#pep3;
		
		for (my $i=0;$i<=$#pep1 ;$i++){
			if($pep1[$i] eq '*'){
				$stop_idx1=$i;
				last;
			}
		}
		for (my $i=0;$i<=$#pep2 ;$i++){
			if($pep2[$i] eq '*'){
				$stop_idx2=$i;
				last;
			}
		}
		for (my $i=0;$i<=$#pep3 ;$i++){
			if($pep3[$i] eq '*'){
				$stop_idx3=$i;
				last;
			}
		}
		if ($stop_idx1>=$stop_idx2 && $stop_idx1>=$stop_idx3){
			$pep=$pep1;
$frame_used=1;

		}
		elsif($stop_idx2>=$stop_idx1 && $stop_idx2>=$stop_idx3){
			$pep=$pep2;
$frame_used=2;
		}
		elsif($stop_idx3>=$stop_idx1 && $stop_idx3>=$stop_idx2){
			$pep=$pep3;
$frame_used=3;
		}
		else{
			print STDERR "ARRRRRRRRRRRRRRRRRRGH!!!!\n";
		}
	}
	else{
		$pep=`transeq -frame $first_frame -auto -stdout -osformat2 plain $nuc_out:$transcript`;
$frame_used=$first_frame;
	}
	$pep=~s/\s//g;
$fasta_remark.=" [strand= $transcripts_strand{$transcript}] [ann_first_frame= $first_frame] [frame_used= $frame_used]";
	print PEP ">$transcript$fasta_remark\n";
	
	for(my $k=0;$k<=length($pep);$k=$k+60){
		print PEP substr($pep,$k,60)."\n";
	}


}
close I; close NUC; close PEP;
#my $pep=`transeq -auto -stdout -osformat2 fasta $nuc_out`;
#print PEP $pep;


#foreach my $transcript (@transcs){
#	my @list=$db->get_feature_by_name('transcript_id' => $vega2gbrowse{$transcript});
#	unless($#list==0){
#		print STDERR "Couldn't find $vega2gbrowse{$transcript} in database\n";
#	}
#	foreach my $feat (@list){
#		my @subfeats= $feat->sub_SeqFeature('CDS');
#		foreach my $subfeat(@subfeats){
#			my $rel_subfeat=$feat->new_from_parent;
#			my $str = $rel_subfeat->gff_string;
#			print "$str\n";
#		}
#	}
#}
