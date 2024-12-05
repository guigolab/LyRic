#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
$|=1;

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";

=head1 NAME

samToPolyA

=head1 SYNOPSIS

A utility to detect poly-adenylated sequencing reads, call on-genome polyA sites and infer the reads' strand based on reads-to-genome alignments in SAM format.

B<Usage example> (on a BAM file):

C<< samtools view $file.bam |samToPolyA.pl --minClipped=20 --minAcontent=0.9  - > ${file}_polyAsites.bed >>


=head2 INPUT

Read-to-genome alignments in SAM format, and the corresponding genome sequence in multifasta.

The script looks for terminal soft-clipped A/T sequences (marked as "S" in the CIGAR string).

=head2 OPTIONS

This script maps polyA sites on the genome based on read mappings in SAM format, and according to the following provided parameters:

=over

=item B<minClipped> (integer) = minimum length of A or T tail required to call a PolyA site.

Default: '10'.

=item B<minAcontent> (float) = required A (or T, if minus strand) content of the tail.

Default: '0.8'.

Note: minAcontent affects both the A tail and the upstream A stretch.

=item B<discardInternallyPrimed> = when enabled, the program will try to avoid outputting false polyA sites arising from internal mis-priming during the cDNA library construction. This option is particularly useful if your cDNA was oligo-dT primed.

Default: disabled.

Requires option B<genomeFasta> to be set.

=item B<minUpMisPrimeAlength> (integer) (ignored if B<discardInternallyPrimed> is not set) = minimum length of genomic A stretch immediately upstream a putative site required to call a false positive (presumably due to internal RT priming), and hence not report the corresponding site in the output.

Default: '10'.

=item B<genomeFasta> (string) (valid only if B<discardInternallyPrimed> is set)= path to multifasta of genome (+ spike-in sequences if applicable), used to extract upstream genomic sequence.

B<Note>: You need write access to the directory containing this file, as the included Bio::DB::Fasta module will create a genomeFasta.index file if it doesn't exist.

=back

=head2 OUTPUT

The script will output BED6 with the following columns:

=over

=item column 1: chromosome

=item column 2: start of polyA site (0-based)

=item column 3: end of polyA site

=item column 4: ID of the read containing a polyA tail

=item column 5: length of the polyA tail on read

=item column 6: genomic strand of the read (see DESCRIPTION below)

=back

=head1 DESCRIPTION

The script will search for read alignment patterns such as:


C<< XXXXXXXXXXXAAAAAAAAAAAAAAA(YYYY) [read] >>

C<< |||||||||||..................... [match] >>

C<< XXXXXXXXXXXZ-------------------- [reference sequence] >>

or

C<< (YYYY)TTTTTTTTTTTTTTTTXXXXXXXXXX [read] >>

C<< ......................|||||||||| [match] >>

C<< ---------------------ZXXXXXXXXXX [reference sequence] >>

Where:

=over

=item C<|> / C<.> = a position mapped / unmapped to the reference, respectively

=item C<X> = the mapped portion of the read or reference sequence

=item C<(Y)> = an optional soft-clipped, non-(A|T)-rich sequence (possibly a sequencing adapter)

=item C<Z> = the position on the reference sequence where the alignment breaks

=item The C<A> / C<T> streches are soft-clipped ('S' in CIGAR nomenclature) in the alignment

=item C<-> = the portion of the reference sequence unaligned to the read

=back

The genomic strand of the read + polyA site is inferred from the mapping of the read, I<i.e.>, reads where a polyA tail was detected at their 3' end are assigned a '+' genomic strand, whereas reads with a polyT tail at their 5' end are deduced to originate from the '-' strand. In that example, the first / second alignment would lead to a called polyA site at position Z on the '+' / '-' strand of the reference sequence, respectively.

=head1 DEPENDENCIES

CPAN: Bio::DB::Fasta

=head1 AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

=cut

my $minSoftClippedSeqLength=10;
my $minAcontent=0.8;
my $minUpMisPrimeAlength=10;
my $discardInternallyPrimed='';
my $genomeFa;
GetOptions ('minClipped=i' => \$minSoftClippedSeqLength,
            'minAcontent:f' => \$minAcontent,
            'minUpMisPrimeAlength=i' => \$minUpMisPrimeAlength,
            'genomeFasta=s' => \$genomeFa,
            'discardInternallyPrimed' => \$discardInternallyPrimed)
  or pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );

unless(defined $minSoftClippedSeqLength && defined $minAcontent){
	pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
}

my $chrdb;

if($discardInternallyPrimed){
	print STDERR "Will try to discard false polyA tails resulting from internal mis-priming.\n";
	if(defined $genomeFa && defined $minUpMisPrimeAlength){
		$chrdb = Bio::DB::Fasta->new("$genomeFa", -reindex => 0);
	}
	else{
		pod2usage( { -message => "Error in command line arguments: --discardInternallyPrimed requires --genomeFasta to be set.\n",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
	}
}

else{
#	print STDERR "Will not try to discard false polyA tails resulting from internal mis-priming.\n";
	if(defined $genomeFa){
		pod2usage( { -message => "Error in command line arguments: --genomeFasta requires --discardInternallyPrimed to be set.\n",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
	}
}

while (<STDIN>){
	my $line=$_;
	chomp;
	next if($_=~/^\@SQ/); #skip sequence header
	my @line=split "\t";
	die "Invalid format (doesn't look like SAM)\n" unless ($#line>9);
	next if($line[5] eq '*'); #skip unmapped reads
	my @cigarNumbers=split (/[A-Z]/,$line[5]);
	my @cigarLetters=split(/\d+/,$line[5]);
	shift(@cigarLetters); #the first element is empty, remove it
	#check if there's any soft-clipped sequence at the end(s) of the read (i.e. unmapped tail)
	if($cigarLetters[0] eq 'S' && $cigarNumbers[0] > $minSoftClippedSeqLength){ #tail is likely at the beginning of the read
	#it means strand is -
		my $tailSeq=substr($line[9], 0, $cigarNumbers[0]);
		#reverse string so that we start scanning the sequence where it drops from the genome
		$tailSeq=reverse($tailSeq);
		$tailSeq=~ tr/ACGTacgt/TGCAtgca/;
		my $strand='-';
		#there might be an exogenous sequence adapter at the end of the sequence.
		my $aTailLength=countAs($tailSeq, $minSoftClippedSeqLength);
		if ($aTailLength>0){
			my $upstreamGenomeAsLength=0;

			if($discardInternallyPrimed){
				my $upstreamStop=$line[3]+$minUpMisPrimeAlength;
				my $upstreamGenomeSeq=$chrdb->seq($line[2], $line[3], $upstreamStop);
				$upstreamGenomeAsLength=0;
				if(defined $upstreamGenomeSeq && length($upstreamGenomeSeq)>0){
					$upstreamGenomeSeq=reverse($upstreamGenomeSeq);
					$upstreamGenomeSeq=~ tr/ACGTacgt/TGCAtgca/;
					$upstreamGenomeAsLength=countAs($upstreamGenomeSeq, $minUpMisPrimeAlength);
				}
				else{
					warn "WARNING: Could not extract reference sequence in $line[2] (is it present in $genomeFa?). Assuming no A stretch upstream of putative site.\n"
				}
			}
			unless ($upstreamGenomeAsLength >= $minUpMisPrimeAlength){
				my $start=$line[3]-1;
				print "$line[2]\t$start\t$line[3]\t$line[0]\t$aTailLength\t$strand\n";
				next; # found, no need to look at the other end of the read
			}
		}
	}

	if($cigarLetters[$#cigarLetters] eq 'S' && $cigarNumbers[$#cigarNumbers] > $minSoftClippedSeqLength){ #tail is likely at the end of the read
	# this is no ELSIF!! need to check second end if first fails
		my $tailSeq=substr($line[9], length($line[9]) - $cigarNumbers[$#cigarNumbers], $cigarNumbers[$#cigarNumbers]);
		my $strand='+';
		#there might be an exogenous sequence adapter at the end of the sequence.
		my $aTailLength=countAs($tailSeq, $minSoftClippedSeqLength);
		if($aTailLength>0){
			my $upstreamGenomeAsLength=0;
				#calculate where the start of the tail is on the genome
			my $genomicLength=0;
			for (my $i=0; $i<=$#cigarLetters;$i++){
				if($cigarLetters[$i] =~ /[MDNXP]/){
					$genomicLength+=$cigarNumbers[$i];
				}
			}
			my $start=$line[3]+$genomicLength-1;
			my $end=$start+1;
			if($discardInternallyPrimed){
				my $upstreamStop=$start;
				my $upstreamStart=$start-$minUpMisPrimeAlength;
				my $upstreamGenomeSeq=$chrdb->seq($line[2], $upstreamStart, $upstreamStop);
				my $upstreamGenomeAsLength=0;
				if(defined $upstreamGenomeSeq && length($upstreamGenomeSeq)>0){
					$upstreamGenomeAsLength=countAs($upstreamGenomeSeq, $minUpMisPrimeAlength);
				}
				else{
					warn "WARNING: Could not extract reference sequence in $line[2] (is it present in $genomeFa?). Assuming no A stretch upstream of putative site.\n"
				}
			}
			unless ($upstreamGenomeAsLength >= $minUpMisPrimeAlength){
				print "$line[2]\t$start\t$end\t$line[0]\t$aTailLength\t$strand\n";
			}
		}
	}
}

sub countAs{
	my $minSeqLength=$_[1];
	my $seq=lc($_[0]);
	my @seq=split("", $seq);
	my $countAs=0;
	my $maxMismatches=$minSeqLength-($minAcontent*$minSeqLength);
	my $SeqLength=0;
	my $countMismatches=0;
	for(my $i=0; $i<=$#seq;$i++){
	 my $nt=$seq[$i];
	 $SeqLength=$i+1;
	 if($countMismatches>$maxMismatches && $SeqLength>$minSeqLength){
		 last;
	 }
	 if($nt eq "a"){
	 	$countAs++;
	 }
	 else{
	 	$countMismatches++;
	 }
	}
	if($SeqLength >= $minSeqLength && $countAs > $minSeqLength-$maxMismatches){
		return $SeqLength;
	}
	else{
		return 0;
	}

}
