#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib "/users/rg/jlagarde/julien_utils/";
use gffToHash;
use overlap;
$|=1;
my $chrAnnotationFile=$ARGV[0];
my $transcriptCoordFile=$ARGV[1];
$Data::Dumper::Sortkeys = 1;

###########################################################################
###########################################################################
####
#### This scripts converts feature coordinates on spliced transcripts to
#### global genomic coordinates.
####
####
##########################
###################### I/O :
##########################
#### Takes 2 GFF files as input, and writes genome GFF on file2 in standard output
####     Arg 1 (file1) : genome annotation (transcripts on genome, chromosome coordinates)
####     Arg 2 (file2) : transcript annotation (features on spliced transcript coordinates)
#### file1 and file2 are joined on the transcript_id GFF attribute (9th field) and seqname (1st field), respectively.
####
####     ###         Input example: file1         ###
####
#### chr10   HAVANA  exon    6622381 6622889 .       +       .       transcript_id "ENST00000445427.1";
#### chr10   HAVANA  exon    6625638 6626261 .       +       .       transcript_id "ENST00000445427.1";
####
####     ###         Input example: file2         ###
####
#### ENST00000445427.1       Primer    Primer    639     665     .       +       .       transcript_id ENST00000445427.1.3p.686;
####     ###         Output example               ###
####
#### chr10     Primer    exon    6625767 6625793 .       +       .       transcript_id "ENST00000445427.1.3p.686";
####
####
##########################
###################### EXAMPLE USAGE:
##########################
####
####     $ transcriptCoords2chromCoords.pl transcripts_on_genome.gff features_on_transcript.gff > features_on_genome.gff
####
###########################################################################
###########################################################################

print STDERR "Reading genome coords...\n";
my %transcriptCoordsGffHash=gffToHash($ARGV[1], 'transcript_id');
print STDERR "Done.\n";
print STDERR "Reading on-transcript coords...\n";
my %chrAnnotationGffHash=gffToHash($ARGV[0], 'transcript_id');
print STDERR "Done.\n";
my %exonsOnTranscript=(); #map exon coordinates to transcript

foreach my $primerId (keys %transcriptCoordsGffHash){
	foreach my $primerRecord (@{$transcriptCoordsGffHash{$primerId}}){
		my $transcriptId=${$primerRecord}[0];
		my $primerTranscriptStart=${$primerRecord}[3];
		my $primerTranscriptEnd=${$primerRecord}[4];
		my $primerTranscriptStrand=${$primerRecord}[6];
		if(exists $chrAnnotationGffHash{$transcriptId}){
			my $transcriptLength=0;
			my $transcriptGenomeStrand='';
			my $transcriptChr='';
			foreach my $gffRecord (@{$chrAnnotationGffHash{$transcriptId}}){
				if(${$gffRecord}[2] eq 'exon'){
					$transcriptLength+=(${$gffRecord}[4]-${$gffRecord}[3])+1;
					$transcriptGenomeStrand=${$gffRecord}[6];
					$transcriptChr=${$gffRecord}[0];
				}
			}
			die "OUT OF COORDINATES: primer $primerId, transcript $transcriptId, transcriptLength $transcriptLength, primerStart $primerTranscriptStart, primerEnd $primerTranscriptEnd\n" if($primerTranscriptStart<1 || $primerTranscriptEnd<1 || $primerTranscriptStart>$transcriptLength || $primerTranscriptEnd>$transcriptLength);
			if(${${$chrAnnotationGffHash{$transcriptId}}[0]}[6] eq '-'){
				#then reverse coordinates on spliced transcript
				my $revPrimerTranscriptStart=($transcriptLength - $primerTranscriptEnd) +1;
				my $revPrimerTranscriptEnd=($transcriptLength - $primerTranscriptStart) +1;
				$primerTranscriptStart=$revPrimerTranscriptStart;
				$primerTranscriptEnd=$revPrimerTranscriptEnd;
			}
			#map exon coordinates to transcript :
 			my $previousExonsLength=0;
			foreach my $gffRecord (sort { $a->[3] <=> $b->[3] } @{$chrAnnotationGffHash{$transcriptId}}){ #sort 	exon records by exon start
				if(${$gffRecord}[2] eq 'exon'){
					my $exonGenomeStart=${$gffRecord}[3];
					my $exonGenomeEnd=${$gffRecord}[4];
					my $exonTranscriptStart=$previousExonsLength+1;
					my $exonTranscriptEnd=$exonTranscriptStart+($exonGenomeEnd-$exonGenomeStart);
					$previousExonsLength=$exonTranscriptEnd;
					@{$exonsOnTranscript{$transcriptId}{$exonGenomeStart}{$exonGenomeEnd}}=($exonTranscriptStart, $exonTranscriptEnd);
				}
			}
			#compare primer coordinates to junctions on transcript coordinates, then deduce chr coordinates
			 my @primerChrStarts=();
			 my @primerChrEnds=();
			 my $primerChrStrand='';
			 foreach my $exonGenomeStart (sort { $a <=> $b } keys %{$exonsOnTranscript{$transcriptId}}){
			 	foreach my $exonGenomeEnd (keys %{$exonsOnTranscript{$transcriptId}{$exonGenomeStart}}){
			 		my $exonTranscriptStart=${$exonsOnTranscript{$transcriptId}{$exonGenomeStart}{$exonGenomeEnd}}[0];
			 		my $exonTranscriptEnd=${$exonsOnTranscript{$transcriptId}{$exonGenomeStart}{$exonGenomeEnd}}[1];
			 		my @overlap=overlap::overlap($primerTranscriptStart,$primerTranscriptEnd,$exonTranscriptStart, $exonTranscriptEnd, $transcriptId, $transcriptId, '+', '+');
			 		if($overlap[0] > 0){ # i.e. they do overlap
			 			if($primerTranscriptStart<$exonTranscriptStart){
			 				push(@primerChrStarts, $exonGenomeStart);
			 			}
			 			elsif($primerTranscriptStart >= $exonTranscriptStart){
			 				push(@primerChrStarts, $exonGenomeStart+($primerTranscriptStart-$exonTranscriptStart))
			 			}
			 			else{
			 				die;
			 			}
			 			if($primerTranscriptEnd <= $exonTranscriptEnd){
			 				push(@primerChrEnds, $exonGenomeStart+($primerTranscriptEnd-$exonTranscriptStart));
			 			}
			 			elsif($primerTranscriptEnd > $exonTranscriptEnd){
			 				push(@primerChrEnds, $exonGenomeEnd);
			 			}
			 			else{
			 				die;
			 			}
			 		}

			 	}
			}
			if($primerTranscriptStrand eq '+' || $primerTranscriptStrand eq '.'){
				if($transcriptGenomeStrand eq '+' || $transcriptGenomeStrand eq '.'){
					$primerChrStrand='+';
				}
				else{
					$primerChrStrand='-';
				}
			}
			else{
				if($transcriptGenomeStrand eq '+' || $transcriptGenomeStrand eq '.'){
					$primerChrStrand='-';
				}
				else{
					$primerChrStrand='+';
				}
			}
			for(my $i=0; $i<=$#primerChrStarts; $i++){
				my @attrs=();
				print "$transcriptChr\t${$primerRecord}[1]\texon\t$primerChrStarts[$i]\t$primerChrEnds[$i]\t.\t$primerChrStrand\t.\t";
				foreach my $key (keys %{${$primerRecord}[8]}){
					if($key eq 'transcript_id'){
					unshift(@attrs, "$key \"${${$primerRecord}[8]}{$key}\";");
					}
					else{
					push(@attrs, "$key \"${${$primerRecord}[8]}{$key}\";");
					}
				}
				print join(" ", @attrs)."\n";
			}

		}
		else{ # i.e. not exists $chrAnnotationGffHash{$transcriptId}
 			warn "transcript_id $transcriptId not found in $ARGV[0]\n."
 		}
	}
}



#############################################################################################
# -spliced_coordinates  34 31      30 27    26  22  21  17    16  12   11  7   6    1
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
# +spliced_coordinates  1  4       5  8     9   13  14  18    19  23   24  28  29   34
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
#                       agtc-------cgat-----atcga---tcccg-----gatct----ataca---tgatct
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
# chr_coordinates       100|       111|     120 124 128 132   138 142  147 151 155  160
#                          103        114

#                                   gat-----atcga---tcccg-----gat
# +3pprimer_spliced_coords          6                           21   +
# +3pprimer_chr_coords              112                         140  +
#                                   cta-----tagct---agggc-----cta
# +5pprimer_spliced_coords          6                           21   -
# +5pprimer_chr_coords              112                         140  -

#                                   gat-----atcga---tcccg-----gat
# -3pprimer_spliced_coords          29                          14   +
# -3pprimer_chr_coords              112                         140  -
# -5pprimer_spliced_coords          29                          14   -
# -5pprimer_chr_coords              112                         140  +
#############################################################################################




