#!/usr/bin/env perl
use strict;
use warnings;
use FindBin; # find present script
use lib "$FindBin::Bin"; # include script's directory (where processJsonToHash.pm is)
use Data::Dumper;
#use lib "/users/rg/jlagarde/julien_utils/";
use gffToHash;
use hashToGff;
#use Storable qw(dclone);

# will add an "exonRank" attribute to each gff record given as input
#exonRank corresponds to the rank of the exon from 5' to 3' of the mRNA (first exon has rank 0), divided by X, where X is the total number of exons of the transcript
#example: exonRank "0";
#one can easily deduce if the exon is 5'-most (exonRank==0) or 3'-most (exonRank==1).WARNING: for monoexonic transcripts exonRank is set to 0 by convention

my $gff_in=$ARGV[0];


print STDERR "Parsing input GFF \n";
my %gffHash=gffToHash($gff_in, 'transcript_id', 0, 'exon');
print STDERR "Done\n";
#print STDERR Dumper \%gffHash;

#my %outGff=dclone(%gffHash); #we'll only edit outGff

#print STDERR Dumper \%gffHash;


foreach my $transcript (keys %gffHash){
	my $firstRec=${$gffHash{$transcript}}[0]; #extract first record of transcript id just to deduce strand
	my $strand=$$firstRec[6];
	my $numberOfExons=$#{$gffHash{$transcript}};
	#print "$transcript\t$numberOfExons\n";
	my $exonRank=0;
	if( $strand eq '+' ){ # sort by increasing start
		foreach my $exonRecord (sort { $a->[3] <=> $b->[3] } @{$gffHash{$transcript}}){
			if($numberOfExons > 0){ #multi-exonic
				${${$exonRecord}[8]}{'exonRank'}=$exonRank/$numberOfExons;
				$exonRank++;
			}
			else{ # monoexonic transcript
				${${$exonRecord}[8]}{'exonRank'}=0;

			}
		}
	}
	elsif($strand eq '-'){ # sort by decreasing start
		foreach my $exonRecord (sort { $b->[3] <=> $a->[3] } @{$gffHash{$transcript}}){ # adjust end coordinates of first exon
			if($numberOfExons > 0){ #multi-exonic
				${${$exonRecord}[8]}{'exonRank'}=$exonRank/$numberOfExons;
				$exonRank++;
			}
			else{ # monoexonic transcript
				${${$exonRecord}[8]}{'exonRank'}=0;
			}

		}
	}
	else{
		die;
	}
}

my @outGff=hashToGff(\%gffHash);
print join("", @outGff);
