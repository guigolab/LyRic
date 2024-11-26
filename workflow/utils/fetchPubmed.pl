#!/usr/bin/env perl

use strict;
use warnings;
use LWP::Simple;
use XML::LibXML;
use Data::Dumper;

my $year=$ARGV[0];
my $search_term=$ARGV[1];

# inspired from https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Basic_Pipelines

# Download PubMed records that are indexed in MeSH for both asthma and
# leukotrienes and were also published in 2009.

my $db = 'pubmed';

my $query = "$search_term"."+AND+$year"."[pdat]";

my $totalPerYearQuery= "$year"."[pdat]";

print STDERR "Query is: $query\n";


#assemble the esearch URL
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#$url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

########## extract total number of paper for year $year

my $url = $base . "egquery.fcgi?&term=$totalPerYearQuery";

#post the esearch URL
my $output = get($url);

#print $output;

my $parser = XML::LibXML->new();
my $respDom = $parser->parse_string($output);

my $total=$respDom->getDocumentElement->findnodes('/Result/eGQueryResult//ResultItem[DbName="pubmed"]/Count')->[0]->to_literal;




########## extract number of papers matching query for year $year


$url = $base . "egquery.fcgi?&term=$query";

#post the esearch URL
$output = get($url);

#print $output;

$parser = XML::LibXML->new();
$respDom = $parser->parse_string($output);

my $count=$respDom->getDocumentElement->findnodes('/Result/eGQueryResult//ResultItem[DbName="pubmed"]/Count')->[0]->to_literal;

my $perHundredThousand=($count/$total)*100000;

print "$year\t$search_term\t$perHundredThousand\n";


sleep 1; # max 3 NCBI queries per second
