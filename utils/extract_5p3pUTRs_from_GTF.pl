#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
$|=1;


# takes a *sorted* GTF as input (stdin), and uses 'CDS' + 'exon' records to build 5pUTR and 3pUTR records
# e.g.: zcat gen10.gtf.gz |sortgff |grep "transcript_type \"protein_coding\";"| extract_5p3pUTRs_from_GTF.pl

my %transcriptAttr=();
my %coding_transcript=();
my %transcriptCdsStarts=();
my %transcriptCdsEnds=();
my %transcriptExonStarts=();
my %transcriptExonEnds=();
my %transcriptStrand=();
my %transcriptChr=();
my %transcriptSrc=();
my %transcript5pUtrStart=();
my %transcript5pUtrEnd=();
my %transcript3pUtrStart=();
my %transcript3pUtrEnd=();

while(<>){
	chomp;
	my @line=split "\t";
	my $line=$_;
	$line[8]=~/transcript_id \"(\S+)\";/;
	my $transcript_id=$1;
	if($line[2] eq 'exon'){
		push (@{$transcriptExonStarts{$transcript_id}}, $line[3]);
		push (@{$transcriptExonEnds{$transcript_id}}, $line[4]);

	}
	elsif($line[2] eq 'CDS'){
		$transcriptAttr{$transcript_id}=$line[8];
		$coding_transcript{$transcript_id}=1;
		$transcriptStrand{$transcript_id}=$line[6];
		$transcriptChr{$transcript_id}=$line[0];
		$transcriptSrc{$transcript_id}=$line[1];
		push (@{$transcriptCdsStarts{$transcript_id}}, $line[3]);
		push (@{$transcriptCdsEnds{$transcript_id}}, $line[4]);
	}
}


foreach my $transcript_id (keys %coding_transcript){
	my @leftUtrStarts;
	my @leftUtrEnds;
	my @rightUtrStarts;
	my @rightUtrEnds;
	# print $transcript_id."\ttranscriptExonStarts\t".join(" , ", @{$transcriptExonStarts{$transcript_id}})."\n";
	# print $transcript_id."\ttranscriptExonEnds\t".join(" , ", @{$transcriptExonEnds{$transcript_id}})."\n";
	# print $transcript_id."\ttranscriptCdsStarts\t".join(" , ", @{$transcriptCdsStarts{$transcript_id}})."\n";
	# print $transcript_id."\ttranscriptCdsEnds\t".join(" , ", @{$transcriptCdsEnds{$transcript_id}})."\n";

	for (my $i=0; $i<=$#{$transcriptExonStarts{$transcript_id}}; $i++){ #loop over exons
		if(${$transcriptExonStarts{$transcript_id}}[$i]<${$transcriptCdsStarts{$transcript_id}}[0]){ #build left UTR
			push(@leftUtrStarts, ${$transcriptExonStarts{$transcript_id}}[$i]);
			if(${$transcriptExonEnds{$transcript_id}}[$i]<${$transcriptCdsStarts{$transcript_id}}[0]){
				push(@leftUtrEnds, ${$transcriptExonEnds{$transcript_id}}[$i]);
			}
			else{
				push(@leftUtrEnds, ${$transcriptCdsStarts{$transcript_id}}[0]-1);
			}
		}
		if(${$transcriptExonEnds{$transcript_id}}[$i]>${$transcriptCdsEnds{$transcript_id}}[$#{$transcriptCdsEnds{$transcript_id}}]){ #build right UTR
			push(@rightUtrEnds, ${$transcriptExonEnds{$transcript_id}}[$i]);
			if(${$transcriptExonStarts{$transcript_id}}[$i]>${$transcriptCdsEnds{$transcript_id}}[$#{$transcriptCdsEnds{$transcript_id}}]){
				push(@rightUtrStarts, ${$transcriptExonStarts{$transcript_id}}[$i]);
			}
			else{
				push(@rightUtrStarts, ${$transcriptCdsEnds{$transcript_id}}[$#{$transcriptCdsEnds{$transcript_id}}]+1);
			}
		}
	}
	# print $transcript_id."\tleftUtrStarts\t".join(" , ", @leftUtrStarts)."\n";
	# print $transcript_id."\tleftUtrEnds\t".join(" , ", @leftUtrEnds)."\n";
	# print $transcript_id."\trightUtrStarts\t".join(" , ", @rightUtrStarts)."\n";
	# print $transcript_id."\trightUtrEnds\t".join(" , ", @rightUtrEnds)."\n";
	my @fivePrimeUtrStarts;
	my @fivePrimeUtrEnds;
	my @threePrimeUtrStarts;
	my @threePrimeUtrEnds;
	if($transcriptStrand{$transcript_id} eq '-'){
		@fivePrimeUtrStarts=@rightUtrStarts;
		@fivePrimeUtrEnds=@rightUtrEnds;
		@threePrimeUtrStarts=@leftUtrStarts;
		@threePrimeUtrEnds=@leftUtrEnds;
	}
	else{
		@fivePrimeUtrStarts=@leftUtrStarts;
		@fivePrimeUtrEnds=@leftUtrEnds;
		@threePrimeUtrStarts=@rightUtrStarts;
		@threePrimeUtrEnds=@rightUtrEnds;
	}

	for (my $i=0; $i<=$#fivePrimeUtrStarts;$i++){
		print "$transcriptChr{$transcript_id}\t$transcriptSrc{$transcript_id}\t5pUTR\t$fivePrimeUtrStarts[$i]\t$fivePrimeUtrEnds[$i]\t.\t$transcriptStrand{$transcript_id}\t.\t$transcriptAttr{$transcript_id}\n";
	}
	for (my $i=0; $i<=$#threePrimeUtrStarts;$i++){
		print "$transcriptChr{$transcript_id}\t$transcriptSrc{$transcript_id}\t3pUTR\t$threePrimeUtrStarts[$i]\t$threePrimeUtrEnds[$i]\t.\t$transcriptStrand{$transcript_id}\t.\t$transcriptAttr{$transcript_id}\n";
	}
}
