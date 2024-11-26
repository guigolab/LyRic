#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use lib "/users/rg/jlagarde/julien_utils/";
use gffToHash;
use hashToGff;
use Storable qw(dclone);

#adds 'transcript' records to 'exon'-only gff

my $transcripts_gff=$ARGV[0];
my %gffHash=gffToHash($transcripts_gff, 'transcript_id', 1, 'exon');

my %transcriptRecordsHash=();

foreach my $tx (keys %gffHash){
	my $firstRec=${$gffHash{$tx}}[0]; #extract first record of transcript id just to deduce strand
	my $strand=$$firstRec[6];
	my $chr=$$firstRec[0];
	my $src=$$firstRec[1];
	my $txStart;
	my $txEnd;
	my $geneId=$$firstRec[8]{'gene_id'};
#	if( $strand eq '+'){
	my @sortedExons=sort { $a->[3] <=> $b->[3] } @{$gffHash{$tx}};
	$txStart=${$sortedExons[0]}[3];
	$txEnd=${$sortedExons[$#sortedExons]}[4];
#	print "$txStart $txEnd\n";
	my @tmpArr = (
            $chr,
            $src,
            'transcript',
            $txStart,
            $txEnd,
            0,
            $strand,
            '.',
            {
              gene_id => $geneId,
              transcript_id => $tx
            }
);
	${$transcriptRecordsHash{$tx}}[0]= \@tmpArr;

}

my @outGff=hashToGff(\%transcriptRecordsHash);
push(@outGff, hashToGff(\%gffHash));
#print Dumper \%transcriptRecordsHash;

print join("", @outGff);