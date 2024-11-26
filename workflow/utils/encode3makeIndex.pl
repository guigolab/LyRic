#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
#use JSON;
use POSIX;
use processJsonToHash;
use lib "$ENV{'ENCODE3_PERLMODS'}";
use lib "/users/rg/jlagarde/julien_utils/";
use indexHashToIndexFile;
use indexFileToHash;
#use processJsonToHash;
#use JSON;
use fetchMetadataFromCollection;
$Data::Dumper::Sortkeys =1;
$Data::Dumper::Deepcopy = 1;
$Data::Dumper::Purity =1;
binmode STDOUT, ":encoding(UTF-8)";
$|=1;
#my $inputList=$ARGV[0]; # path to file containing the list of accessions or uuid's to search
# one per line. Quotes optional.
# The script will then look into entire objects collection for the object with accession=XXX
my $propertiesToOutput=$ARGV[0]; # path to file
my $outputFormat=$ARGV[1];
my $objectType=$ARGV[2];

unless ($#ARGV==2){
	die "Wrong dumber of args. Should be 3. Exiting.\n"
}
unless ($ARGV[1] eq 'index' || $ARGV[1] eq 'tsv'){
	die "arg2 should be 'index' or 'tsv'. Exiting.\n"
}

#my $encodeFilesBasedir=$ENV{'ENCODE3_FILE_REPOSITORY'}."/www.encodeproject.org";
my $encodeFilesBasedir="https://www.encodeproject.org";

#open KEYS, "$inputList" or die "$inputList: $!";
my $collectionJsonFile=$ENV{'ENCODE3_FULL_OBJECTS_COLLECTION'};
#print STDERR "Converting $collectionJsonFile to hash...";
my $fullCollection = processJsonToHash($collectionJsonFile);
#print STDERR " Done.\n";
#print STDERR Dumper $fullCollection;

my %keys=();

foreach my $id (keys %{$fullCollection}){
	 foreach my $i ( @{$$fullCollection{$id}{'@type'}}){
	 	if ( $i eq $objectType){
	 		if(exists $$fullCollection{$id}{'accession'}){
	 			$keys{$$fullCollection{$id}{'accession'}}=1;
	 			last;
	 		}
	 	}
	 }
}

#print STDERR Dumper \%keys;

#print STDERR "fullColl ".Dumper %{$fullCollection}."\nend fullColl\n";
#%$fullCollection=(); #this is useless, as it doesn't free up any RAM
#print STDERR "fullColl ".Dumper %{$fullCollection}."\nend fullColl\n";


open PROP, "$propertiesToOutput" or die "$propertiesToOutput: $!";
my %properties=();

while(<PROP>){
	next unless $_=~/\S+/;
	next if ($_=~/^#/);
	chomp;
	my $line=$_;
	$properties{$line}=1;
}
close PROP;

my $metadataHash=fetchMetadataFromCollection(\%keys, \%properties);
#print STDERR Dumper $metadataHash;

#{
#no warnings 'uninitialized';
if($outputFormat eq 'index'){
	my %metadataHashIndexedByFilePath=();
	foreach my $key (keys %{$metadataHash}){
		my $newKey=$encodeFilesBasedir.${$metadataHash}{$key}{'href'};
		$metadataHashIndexedByFilePath{$newKey}=\%{${$metadataHash}{$key}};
	}
	#print STDERR Dumper \%metadataHashIndexedByFilePath;
	my $indexFile=indexHashToIndexFile(\%metadataHashIndexedByFilePath);
	print $indexFile;
}

elsif ($outputFormat eq 'tsv'){
	my $fileDate=POSIX::strftime(
             "%d-%b-%Y",
             localtime(
                 ( stat "$collectionJsonFile" )[9]
                 )
             );

	my @properties=();
	foreach my $prop (sort keys %properties){
		push (@properties, $prop);
	}
	print "#".join("\t", @properties)."\t[ jsonDownloaded: $fileDate ]\n";
	foreach my $key (sort keys %{$metadataHash}){
		my @values=();
		foreach my $prop (@properties){
			my $val='';
			if (defined ($$metadataHash{$key}{$prop})){
				$val=$$metadataHash{$key}{$prop}
			}
			push(@values, $val);
		}
		print join ("\t", @values)."\n";
	}

}

else{
	die;
}
#}
