#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;

# maps AAs to genome sequence from CDS records in input gff files
#- this is based on hg17
#- format is self-explanatory i think. At exon boundaries it will look like:

#AP000295.6-002  72      S       chr21   +       33539248,33540893       33539249,33540893

#- the input used was this file: ftp://genome.imim.es/pub/other/gencode/data/havana-encode/version02.2_14oct05/44regions/44regions_coding_CHR_coord.gtf.
#the mapping were made using only CDS coordinates and the frame info of the first ('first' in the 5'->3' direction of the transcript) CDS exon of each transcript. As I mentioned in http://genome.imim.es/biosapiens/gencode/dataset/v2.2.problems.html the frame info is often erroneous in this release of gencode, which means you'll come accross quite a few in-frame stops in this file!
#Also, please bear in mind that in order to overcome this problem this frame info was NOT used when i generated the AA sequences available at http://genome.imim.es/biosapiens/gencode/dataset/v2.2.test.html (see ftp://genome.imim.es/pub/other/gencode/data/seqs/gencode_analysis/oct2005/CODING_gencode_annotations_44regions_hg17_CHR_coord_oct2005.gtf_CDS.readme)






open IN, "$ARGV[0]" or die $!;

my %transcripts_cdsstarts;
my %transcripts_cdsends;
my %frame;
my %strand;
my %cds_exon_id;
my %chr;

while(<IN>){
	if($_=~/^(\S+)\t(\S+)\tCDS\t(\d+)\t(\d+)\t\S+\t(\S+)\t(\S+)\t.*transcript_id "(\S+)?";.*exon_id "(\S+)?";.*\n/){
		
		push(@{$transcripts_cdsstarts{$7}},$3);
		push(@{$transcripts_cdsends{$7}},$4);
		${$frame{$7}}{$3}{$4}=$6;
		${$cds_exon_id{$7}}{$3}{$4}=$8;
		$strand{$7}=$5;
		$chr{$7}=$1;
	}
}
print STDOUT "#transcript_id\tAA_position\tAA\tcodon\tchr\tstrand\tstart(s)\tstop(s)\n";
foreach my $transcript (keys %transcripts_cdsstarts){
	my $chr_file=$chr{$transcript}.".fa";
	if($strand{$transcript} eq '+'){
		@{$transcripts_cdsends{$transcript}}=sort { $a <=> $b } @{$transcripts_cdsends{$transcript}};
		@{$transcripts_cdsstarts{$transcript}}=sort { $a <=> $b } @{$transcripts_cdsstarts{$transcript}};
	}
	else {
		@{$transcripts_cdsends{$transcript}}=sort { $b <=> $a } @{$transcripts_cdsends{$transcript}};
		@{$transcripts_cdsstarts{$transcript}}=sort { $b <=> $a } @{$transcripts_cdsstarts{$transcript}};
	}
	#print STDERR "\n".$transcript;
	my $aapos=0;
	for (my $i=0; $i <=$#{$transcripts_cdsstarts{$transcript}};$i++){
		#print STDERR "\n\tCDS ${$transcripts_cdsstarts{$transcript}}[$i] ${$transcripts_cdsends{$transcript}}[$i]";
		
		if($strand{$transcript} eq '+'){
			for(my $j=${$transcripts_cdsstarts{$transcript}}[$i]+$frame{$transcript}{${$transcripts_cdsstarts{$transcript}}[$i]}{${$transcripts_cdsends{$transcript}}[$i]}; $j<${$transcripts_cdsends{$transcript}}[$i];$j=$j+3){
				$aapos++;
				my @aastarts=();
				my @aastops=();
				my $codon='';
				my $aa='';
				#	print STDOUT "\n\t\t$j, ${$transcripts_cdsends{$transcript}}[$i]\n";
				if($j+2>${$transcripts_cdsends{$transcript}}[$i]){
					if (exists ${$transcripts_cdsstarts{$transcript}}[$i+1]){
						@aastarts=($j, ${$transcripts_cdsstarts{$transcript}}[$i+1]);
						@aastops=(${$transcripts_cdsends{$transcript}}[$i], ${$transcripts_cdsstarts{$transcript}}[$i+1]+ $j+1-${$transcripts_cdsends{$transcript}}[$i]);
						
						if($j+2-${$transcripts_cdsends{$transcript}}[$i] != $frame{$transcript}{${$transcripts_cdsstarts{$transcript}}[$i+1]}{${$transcripts_cdsends{$transcript}}[$i+1]}){
							my $calc=$j+2-${$transcripts_cdsends{$transcript}}[$i];
							print STDERR "# Inconsistent frame (ann. $frame{$transcript}{${$transcripts_cdsstarts{$transcript}}[$i+1]}{${$transcripts_cdsends{$transcript}}[$i+1]}, calc. $calc  for $transcript CDS ${$transcripts_cdsstarts{$transcript}}[$i+1] -> ${$transcripts_cdsends{$transcript}}[$i+1])\n";
						}
					}
				}
				else{
					@aastarts=($j);
					@aastops=($j+2);
				}
				my $listaastarts= join ",", @aastarts;
				my $listaastops=join ",", @aastops;
				for(my $k=0; $k<=$#aastarts;$k++){
					$codon.=`chr_subseq /seq/genomes/H.sapiens/golden_path_200405/chromFa/$chr_file $aastarts[$k] $aastops[$k]`;
					chomp $codon;
				}
				
				$aa=`echo $codon | transeq -frame 1 -auto -filter -osformat2 plain`;
				chomp $aa;
				$codon=uc($codon);
				my $line="$transcript\t".$aapos."\t$aa\t$codon\t$chr{$transcript}\t$strand{$transcript}\t$listaastarts\t$listaastops\n";
				if($line=~/^\S+\t\d+\t[A-Z\*]\t[ATGC]{3}\t\S+\t\S\t\S+\t\S+\n/){
					print STDOUT $line;
				}
			}
			
			
		}
		else{ #minus strand
			for(my $j=${$transcripts_cdsends{$transcript}}[$i]-$frame{$transcript}{${$transcripts_cdsstarts{$transcript}}[$i]}{${$transcripts_cdsends{$transcript}}[$i]}; $j>${$transcripts_cdsstarts{$transcript}}[$i];$j=$j-3){
				$aapos++;
				my @aastarts=();
				my @aastops=();
				my $codon='';
				my $aa='';
				if($j-2<${$transcripts_cdsstarts{$transcript}}[$i]){
					if (exists ${$transcripts_cdsends{$transcript}}[$i+1]){
						@aastarts=(${$transcripts_cdsends{$transcript}}[$i+1]- (${$transcripts_cdsstarts{$transcript}}[$i]-$j+1), ${$transcripts_cdsstarts{$transcript}}[$i]);
						@aastops=(${$transcripts_cdsends{$transcript}}[$i+1] , $j);
					}
				}
				else{
					@aastarts=($j-2);
					@aastops=($j);
				}
				my $listaastarts= join ",", @aastarts;
				my $listaastops=join ",", @aastops;
				for(my $k=0; $k<=$#aastarts;$k++){
					$codon.=`chr_subseq /seq/genomes/H.sapiens/golden_path_200405/chromFa/$chr_file $aastarts[$k] $aastops[$k]`;
					chomp $codon;
				}
				$codon=`echo $codon | revseq -auto -filter -osformat2 plain`;
				chomp $codon;
				$aa=`echo $codon | transeq -frame 1 -auto -filter -osformat2 plain`;
				chomp $aa;
				$codon=uc($codon);
				my $line="$transcript\t".$aapos."\t$aa\t$codon\t$chr{$transcript}\t$strand{$transcript}\t$listaastarts\t$listaastops\n";
				if($line=~/^\S+\t\d+\t[A-Z\*]\t[ATGC]{3}\t\S+\t\S\t\S+\t\S+\n/){
					print STDOUT $line;
				}
				
			}
		}
	}
}
