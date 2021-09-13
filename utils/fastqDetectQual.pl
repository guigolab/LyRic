#!/usr/bin/perl -w

# http://wiki.bits.vib.be/index.php/Identify_the_Phred_scale_of_quality_scores_used_in_fastQ

use strict;
use File::Basename;
use List::MoreUtils qw( minmax );

# fastq_detect.pl fastq.file sample-size
# detect fastQ format from quality scores in fastQ input file
# Version 3
#
# Stephane Plaisance - VIB-BITS - July-04-2012
# Joachim Jacob - Aug-02-2012 - joachim.jacob@gmail.com
# - changed the maximum value of Sanger to 73
# - changed reading the file with a file handle
#   (was a file handle !! supporting several archive formats. SP)
# - changed the diagnosing algoritm
# Stephane Plaisance - VIB-BITS - April-08-2013
# - merged both versions and corrected flaw in min/max
# thanks to Sergey Mitrfanov for perl reformatting

#####################################################################
# diagnose
#   SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
#   ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
#   ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
#   .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
#   LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
#   !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
#   |                         |    |        |                              |                     |
#  33                        59   64       73                            104                   126
# S 0........................26...31.......40
# X                          -5....0........9.............................40
# I                                0........9.............................40
# J                                   3.....9.............................40
# L 0.2......................26...31........41
#
#  S - Sanger        Phred+33,  raw reads typically (0, 40)
#  X - Solexa        Solexa+64, raw reads typically (-5, 40)
#  I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
#  J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40) with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
#  L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
#####################################################################

my $script = basename($0);

@ARGV gt 0 or die "usage: $script <fastq file> <opt:sample-size (100)>\n";
my ($inputfile, $limit) = @ARGV;
if (! defined $limit) { $limit = 100}; # check first 100 records

my $cnt=0;
my ($min, $max); # global min and max values

print STDERR "\n## Analysing ".$limit." records from $inputfile ... \n";
my $z = ReadFile ($inputfile) || die "Error: cannot read from variant file $inputfile: $!\n";

## parse
while (my $id = <$z>) {
	$id =~ m/^@/ || die "expected @ not found in line 1!\n";
	my $seq = <$z>;
	my $sep = <$z>;
	$sep =~ m/^\+/ || die "expected + not found in line 3!\n";
	my $qual = <$z>;
	chomp($qual);
	$cnt++;
	$cnt>=$limit && last;

	# char to ascii
	my @chars = split("", $qual);
	my @nums = sort { $a <=> $b } (map { unpack("C*", $_ )} @chars);

	if ($cnt==1) {
		($min, $max) = minmax @nums;
	} else {
		my ($lmin, $lmax) = minmax @nums; # local values for this read
		$lmin<$min ? $min=$lmin : $min=$min;
		$lmax>$max ? $max=$lmax : $max=$max;
	}
}

undef $z;

## diagnose
my %diag=(
			'Sanger'		=> '.',
			'Solexa'		=> '.',
			'Illumina 1.3+'	=> '.',
			'Illumina 1.5+'	=> '.',
			'Illumina 1.8+'	=> '.',
			);

my %comment=(
			'Sanger'		=> 'Phred+33,  Q[33; 73],  (0, 40)',
			'Solexa'		=> 'Solexa+64, Q[59; 104], (-5, 40)',
			'Illumina 1.3+'	=> 'Phred+64,  Q[64; 104], (0, 40)',
			'Illumina 1.5+'	=> 'Phred+64,  Q[66; 104], (3, 40), with 0=N/A, 1=N/A, 2=Read Segment Quality Control Indicator',
			'Illumina 1.8+'	=> 'Phred+33,  Q[33; 74],  (0, 41)',
			);

if ($min<33 || $max>104) { die "Quality values corrupt. found [$min; $max] where [33; 104] was expected\n"; }
if ($min>=33 && $max<=73)  {$diag{'Sanger'}='x';}
if ($min>=59 && $max<=104) {$diag{'Solexa'}='x';}
if ($min>=64 && $max<=104) {$diag{'Illumina 1.3+'}='x';}
if ($min>=66 && $max<=104) {$diag{'Illumina 1.5+'}='x';}
if ($min>=33 && $max<=74)  {$diag{'Illumina 1.8+'}='x';}

## report
print STDERR "# sampled raw quality values are in the range of [".$min."; ".$max."]\n";
print STDERR "# format(s) marked below with 'x' agree with this range\n";

foreach my $format (sort keys %diag) {
	print STDERR sprintf("  %-13s : %2s  [%-30s] \n", $format, $diag{$format}, $comment{$format});
}


##############
#### Subs ####

# reads from uncompressed, gzipped and bgzip fastQ files
sub ReadFile {
	my $infile = shift;
	my $FH;
	if ($infile =~ /.bz2$/) {
		open ($FH, "bzcat $infile |") or die ("$!: can't open file $infile");
	} elsif ($infile =~ /.gz$/) {
		open ($FH, "zcat $infile |") or die ("$!: can't open file $infile");
	} elsif ($infile =~ /.fq|.fastq|.txt$/) {
		open ($FH, "cat $infile |") or die ("$!: can't open file $infile");
	} else {
		die ("$!: do not recognise file type $infile");
	}
	return $FH;
}

