#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use File::Basename;
use Bio::DB::Fasta;
use Getopt::Long;
use Text::LevenshteinXS qw(distance);
#use Text::Levenshtein::XS qw/distance/;
$|=1;
####### output:
##
## *.introns.gff
#  all introns records present in input.

#   Strand (in field 7) is the strand of the transcript, as inferred by the majority of SJ motifs of the transcript it belongs to (is '.' if no canonical junction is present in whole transcript, or in case of tie over all transcript's introns )

# "seq" attribute (9th field) contains the SJ dinucleotide motif (according to strand of transcript)

# "extendedDonSeq" and "extendedAccSeq" attributes (9th field) contains +/- $extendedSeqLength nucleotides around donor/acceptor site. This is used to detect potential RT template switching as opposed to genuine splicing. Intron boundary is marked with a "/" in sequence.
#    example of likely RT template switching "substrate" (duplicate seq in uppercase, SJ marked with "/":
#                donor=taccCAGCCC/Caggcaccag
#                acceptor= aggagCAGCCCC/cagaactcc

# extendedSeqEditDist attribute (9th field) = Levenshtein edit distance between extendedDonSeq and extendedAccSeq. The higher, the less likely an RT template switching event to occur.
#
# intronDupScore is extendedSeqEditDist/$extendedSeqLength (the higher, the better)
##  *.transcripts.tsv
# strand of each transcript present in input
# field1 -> transcript_id
# field2 -> strand, as inferred from SJ motifs (GT/AG vs. CT/AC) of all its introns
# field3 -> strand, as provided in the input (meaningless if cDNA alignments)
# field4 -> numebr of introns in transcripts
# field5 -> number of canonical GT/AG introns in trancript
# field6 -> % canonical introns in transcript

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to

my $extendedSeqLength=6;
my $genomeFa=$ARGV[1];
my $outprefix=basename($ARGV[0]);
my $chrdb = Bio::DB::Fasta->new("$genomeFa", -reindex => 0);
if(defined($ARGV[2]) && $ARGV[2] ne ''){
	$outprefix=$ARGV[2];
}

GetOptions ('extSeqLength=i' => \$extendedSeqLength,
            )
or pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );


my $outIntrons=$outprefix.".introns.gff";
my $outTranscripts=$outprefix.".transcripts.tsv";

open GFF, "$ARGV[0]" or die $!;

my $count_introns=0;
my %transcripts=();
my %introns=();
print STDERR "Reading GFF in and extracting SJ sequences...\n";
while(<GFF>){
#	print "#".$_;
	chomp;
	next if $_=~/^#/;
	my @line=split "\t";
	next unless ($line[2] eq 'intron');
	$count_introns++;
	$line[8]=~/transcript_id \"(\S+)\";/;
	my $transcript_id=$1;
#	print STDERR $transcript_id."\n";
	#my $chr_file=$line[0].".fa";
	my $start=$line[3]-$extendedSeqLength;
	my $end=$line[4]+$extendedSeqLength;
	die "start ($start) is >= stop ($end). Offending line: ".join("\t", @line)."\n" if ($start >= $end);
	#my $seq=`~jlagarde/julien_utils/chr_subseq $chromdir/$chr_file $start $end`;
	#print STDERR "$line[0] $start $end\n";
	my $intronId="$line[0].$line[3].$line[4].$line[6]";
	push(@{$transcripts{$transcript_id}}, $intronId);

	if ($count_introns%10000 == 0){
		print STDERR "\tProcessed $count_introns intron records\n";
	}

	next if (exists $introns{$intronId});
	my $seq=uc($chrdb->seq("$line[0]", $start => $end));
	die "Could not extract sequence for ".join("\t", @line)."\n" unless (defined($seq) && $seq ne '');
	#chomp($seq);
	#print $seq."\n";
	my $leftMotif=substr($seq,$extendedSeqLength,2);
	my $rightMotif=substr($seq, length($seq)-$extendedSeqLength-2,2);
	my $leftMotifExtended=substr($seq,0,$extendedSeqLength*2);
	my $rightMotifExtended=substr($seq, length($seq)-$extendedSeqLength*2);
#	print "leftE $leftMotifExtended\nrightE $rightMotifExtended\n";
	#print "$leftMotif $rightMotif\n";
	my $fivep;
	my $threep;
	my $strand;
	my $fivepExtended;
	my $threepExtended;
	if(($leftMotif eq 'GT' || $leftMotif eq 'GC') && $rightMotif eq 'AG'){ # + strand
		$fivep=$leftMotif;
		$threep=$rightMotif;
		$fivepExtended=$leftMotifExtended;
		$threepExtended=$rightMotifExtended;
		$strand='+';
	}
	elsif($leftMotif eq 'CT' && ($rightMotif eq 'AC' || $rightMotif eq 'GC')){ # - strand
		$threep=reverse($leftMotif);
		$fivep=reverse($rightMotif);
		$fivep=~ tr/ACGTacgt/TGCAtgca/;
		$threep=~ tr/ACGTacgt/TGCAtgca/;
		$threepExtended=reverse($leftMotifExtended);
		$fivepExtended=reverse($rightMotifExtended);
		$fivepExtended=~ tr/ACGTacgt/TGCAtgca/;
		$threepExtended=~ tr/ACGTacgt/TGCAtgca/;
		$strand='-';
	}
	else{
		$threep=$rightMotif;
		$fivep=$leftMotif;
		$threepExtended=$rightMotifExtended;
		$fivepExtended=$leftMotifExtended;

		$strand='.';
	}
#	print "fivepE $fivepExtended\nthreepE $threepExtended\n";

	#print "$fivep $threep $strand\n";
#	my $intronLevenshtein=distance($fivepExtended,$threepExtended);
#	my $intronSeqSimilarityScore=$intronLevenshtein/($extendedSeqLength*2);
	my @intron=($line[0], $line[3], $line[4], $line[6], $fivep, $threep, $strand, $fivepExtended, $threepExtended); #, $intronLevenshtein, $intronSeqSimilarityScore);
	$introns{$intronId} = [@intron];
#	print "$transcript_id: ".join (" ", @intron)."\n";

}
print STDERR "DONE reading GFF in and extracting SJ sequences.\n";


foreach my $intronId (keys %introns){
	my $extendedExonDonSeq=substr($introns{$intronId}[7], 0, $extendedSeqLength);
	my $extendedIntronAccSeq=substr($introns{$intronId}[8], 0, $extendedSeqLength);
	my $extendedIntronDonSeq=substr($introns{$intronId}[7], $extendedSeqLength);
	my $extendedExonAccSeq=substr($introns{$intronId}[8],$extendedSeqLength);

	my $intronExonLevenshtein=distance($extendedIntronDonSeq, $extendedExonAccSeq);
	my $exonIntronLevenshtein=distance($extendedExonDonSeq, $extendedIntronAccSeq);
	my $intronLevenshtein;
	if ($intronExonLevenshtein >= $exonIntronLevenshtein){
		$intronLevenshtein=$exonIntronLevenshtein;
	}
	else{
		$intronLevenshtein=$intronExonLevenshtein;
	}

	my $intronSeqSimilarityScore=$intronLevenshtein/$extendedSeqLength;
	push(@{$introns{$intronId}}, $intronLevenshtein, $intronSeqSimilarityScore);

}
#print Dumper \%introns;



#my %transcript_id_to_AnnStrand=();
#my %transcript_id_to_InfStrand=();
open TX, ">$outTranscripts" or die $!;
open IN, ">$outIntrons" or die $!;
print TX "#transcript_id\tinferredStrand\tannotatedStrandInInput\tnumberOfSpliceJunctions\tnumberOfCanonicalSpliceJunctions\tproportionOfCanonicalSpliceJunctions\tproportionOfHighConfidenceSpliceJunctions\n";
foreach my $transcript_id (keys %transcripts){
	my $proportionCanonical=0;
	my $proportionNonSuspicious=0;
	my $totalIntronCount=0;
	my $canonicalIntronCount=0;
	my $nonSuspiciousIntronCount=0;

	my @infIntronChainStrand;
	my @annIntronChainStrand;
	my $transcript_inferred_strand;
	foreach my $intronId (@{$transcripts{$transcript_id}}){
		push(@infIntronChainStrand, $introns{$intronId}[6]);
		push(@annIntronChainStrand, $introns{$intronId}[3]);
	}
	my $annStrand=$annIntronChainStrand[0];
	my $countPlus=0;
	my $countMinus=0;
	my $countDot=0;
	foreach my $intronStrand (@infIntronChainStrand){
		if($intronStrand eq '+'){
			$countPlus++;
		}
		elsif($intronStrand eq '-'){
			$countMinus++;
		}
		elsif($intronStrand eq '.'){
			$countDot++;
		}
		else{
			die;
		}
	}
	#print "$countPlus + , $countMinus -, $countDot .\n";
	if($countPlus > $countMinus ){
		$transcript_inferred_strand='+';
	}
	elsif($countPlus < $countMinus ){
		$transcript_inferred_strand='-';
	}
	else{
		$transcript_inferred_strand='.';
	}
	my @transcript_introns=();
	# re-strand intron motifs according to overall transcript strand
	foreach my $intronId (@{$transcripts{$transcript_id}}){
		$totalIntronCount++;
		my $threep;
		my $fivep;
		my $threepExtended;
		my $fivepExtended;
		if ($introns{$intronId}[6] ne $transcript_inferred_strand){
			#$introns{$intronId}[6]=$transcript_id_to_InfStrand{$transcript_id};
			$threep=reverse($introns{$intronId}[4]);
			$fivep=reverse($introns{$intronId}[5]);
			$fivep=~ tr/ACGTacgt/TGCAtgca/;
			$threep=~ tr/ACGTacgt/TGCAtgca/;
			#$introns{$intronId}[4]=$fivep;
			#$introns{$intronId}[5]=$threep;

			$threepExtended=reverse($introns{$intronId}[7]);
			$fivepExtended=reverse($introns{$intronId}[8]);
			$fivepExtended=~ tr/ACGTacgt/TGCAtgca/;
			$threepExtended=~ tr/ACGTacgt/TGCAtgca/;
			#$introns{$intronId}[7]=$fivepExtended;
			#$introns{$intronId}[8]=$threepExtended;
		}
		else{
			$threep=$introns{$intronId}[5];
			$fivep=$introns{$intronId}[4];
			$threepExtended=$introns{$intronId}[8];
			$fivepExtended=$introns{$intronId}[7];
		}
		my @intron=($introns{$intronId}[0], $introns{$intronId}[1], $introns{$intronId}[2], $transcript_inferred_strand, $fivep, $threep, $fivepExtended, $threepExtended, $introns{$intronId}[9], $introns{$intronId}[10]);

		if(($intron[4] eq 'GT' || $intron[4] eq 'GC') && $intron[5] eq 'AG'){
			$canonicalIntronCount++;
		}

		if($intron[9]>0.1){

			$nonSuspiciousIntronCount++;
		}

		push(@transcript_introns, [@intron]);
		if($totalIntronCount == 0){
			$proportionCanonical=0;
			$proportionNonSuspicious=0;
		}
		else{
			$proportionCanonical=$canonicalIntronCount/$totalIntronCount;
			$proportionNonSuspicious=$nonSuspiciousIntronCount/$totalIntronCount;
		}

	}

#	print "$transcript_id ". Dumper \@transcript_introns;
	print TX "$transcript_id\t$transcript_inferred_strand\t$annStrand\t$totalIntronCount\t$canonicalIntronCount\t$proportionCanonical\t$proportionNonSuspicious\n";

	foreach my $intron (@transcript_introns){
		my $extendedDonSeq=substr($$intron[6], 0, $extendedSeqLength)."/".substr($$intron[6], $extendedSeqLength);
		my $extendedAccSeq=substr($$intron[7], 0, $extendedSeqLength)."/".substr($$intron[7], $extendedSeqLength);

		print IN "$$intron[0]\tinfIntron\tintron\t$$intron[1]\t$$intron[2]\t.\t$transcript_inferred_strand\t.\ttranscript_id \"$transcript_id\"; gene_id \"$transcript_id\"; intron_id \"$$intron[0]_$$intron[1]_$$intron[2]_$transcript_inferred_strand\"; seq \"$$intron[4]/$$intron[5]\"; extendedDonSeq \"$extendedDonSeq\"; extendedAccSeq \"$extendedAccSeq\"; extendedSeqEditDist \"$$intron[8]\"; intronDupScore \"$$intron[9]\"; \n";
	}




}



print STDERR "(Found $count_introns intron records in file). Output is in:\n$outTranscripts\nand\n$outIntrons\n";
