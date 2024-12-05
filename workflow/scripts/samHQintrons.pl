#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage;
$|=1;

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";

# =head1 AUTHOR

# Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

# =cut



my $meanQualCutoff=0;
my $window=3; #+/- $window of exonic bases will be extracted around each intron

GetOptions ('slop=i' => \$window,
            'minQual=i' => \$meanQualCutoff)
  or pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );

unless(defined $window && defined $meanQualCutoff && $window >=0 ){
	pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
}


REC: while (<>){
	my $line=$_;
	chomp;
	if($_=~/^\@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_=~/^\@CO\t.*/){ #skip sequence header
		print "$line";
		next REC;
	}
	my @line=split "\t";
	die "Invalid format (doesn't look like SAM)\n" unless ($#line>9);
	next if($line[5] eq '*'); #skip unmapped reads
	#print "$line[0]\n";
	my @cigarNumbers=split (/[A-Z]/,$line[5]);
	my @cigarLetters=split(/\d+/,$line[5]);
	shift(@cigarLetters); #the first element is empty, remove it
	my @qualsAscii=split("", $line[10]);
	my @seq=split("", $line[9]);
	my @quals=();
	foreach my $qual (@qualsAscii){
		push(@quals, ord($qual) -33 ); # see https://en.wikipedia.org/wiki/FASTQ_format#Quality
	}
	@qualsAscii=();
	#print @cigarLetters."\n";
	#print @cigarNumbers."\n";
	#print join (" ", @cigarLetters)."\n";
	#print join (" ", @cigarNumbers)."\n";
	#print join (" ", @quals)."\n";

	my $position=0;
	for (my $i=0; $i<=$#cigarLetters; $i++){
		if($cigarLetters[$i] eq 'M' || $cigarLetters[$i] eq 'I' || $cigarLetters[$i] eq 'S' || $cigarLetters[$i] eq '=' || $cigarLetters[$i] eq 'X' ){
			$position+=$cigarNumbers[$i];
		}
		elsif($cigarLetters[$i] eq 'N'){ #intron!
			my $windowMin=$position-$window;
			my $windowMax=$position+$window;
			$windowMin=0 if $windowMin<0;
			$windowMax=$#seq if $windowMax > $#seq;
			#print join ("", @seq[$windowMin..$position-1])." / ". join("",@seq[$position..$windowMax-1])."\n";
			my $sumQualScores=0;
			my $countBases=0;
			foreach my $qual (@quals[$windowMin..$position-1]) {
				$countBases++;
				$sumQualScores+=$qual;
			}
			foreach my $qual (@quals[$position..$windowMax-1]) {
				$countBases++;
				$sumQualScores+=$qual;
			}
			my $meanQual=$sumQualScores/$countBases;
			#print "$meanQual ($sumQualScores / $countBases)\n";
			if ($meanQual<$meanQualCutoff){
				next REC;
			}
		}
	}
	print "$line"; #also print all monoexonic reads
}
