#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use DBI;
use POSIX;
$|=1;
# input is the oputput of e.g. intersectBed -wb -a stdin -b gen10.long.gene.labels.unstranded.gtf:
# chr1    14590   14666   MICHAELJACKSON_0007:6:55:9788:4378#0/1  255     -       chr1    HAVANA  gene    14363   29806   .       -       .       gene_id "ENSG00000227232.3"; transcript_id "ENSG00000227232.3"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";


while (<>){
	#print;
	chomp;
	my $line=$_;
	my @line=split "\t";
	my $midRead=ceil($line[2]-(($line[2]-$line[1])/2));
#	my $midRead=$line[2]-(($line[2]-$line[1])/2);
	my $geneStrand=$line[12];
	my $geneStart=$line[9];
	my $geneEnd=$line[10];
	my $readDistFrom5prime=undef;

	if($geneStrand eq '-'){
		$readDistFrom5prime=$geneEnd-$midRead;
	}
	else{
		$readDistFrom5prime=$midRead-$geneStart;
	}
	#if($readDistFrom5prime<0){
	#	$readDistFrom5prime=0;
	#}
	#elsif($readDistFrom5prime>($geneEnd-$geneStart)){
	#	$readDistFrom5prime=$geneEnd-$geneStart;
	#}
	my $readNormDistFrom5prime=$readDistFrom5prime/($geneEnd-$geneStart);
	print "$line readNormDistTo5prime \"$readNormDistFrom5prime\";\n";
}
