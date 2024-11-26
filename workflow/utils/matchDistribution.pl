#!/usr/bin/env perl
use Getopt::Long;
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use List::Util 'shuffle';
$|=1;
$Data::Dumper::Sortkeys =1;

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";

=head1 NAME

matchDistribution

=head1 SYNOPSIS

Given distinct B<"subject" (S)> and a B<"target" (T)> distributions, this script attempts to mimic T's density (I<i.e.>, its shape) by randomly sampling from bins in S's population.

B<Usage>: C<< matchDistribution.pl <OPTIONS> <arg1> <arg2> <arg3> >>


=head2 ARGUMENTS/INPUT

=over

=item B<arg1>: Path to file containing T's values (1 column, 1 value per line).

=item B<arg2>: Number of bins to split the distributions into.

=item B<arg3>: Path to tab-separated file containing S's identifiers and values (column 1: unique identifier; column 2: corresponding value).

=back

=head2 OPTIONS

=over

=item B<transform> (string)
= bin transform-transformed values in both distributions. Output values will be the original, non-transformed ones, though.

Possible values: 'log10' only.

Note: binning into log10-transformed is highly recommended e.g. for matching FPKM/RPKM distributions.

=item B<verbose>
= make STDERR more verbose

=back

=head2 OUTPUT

To STDOUT.

The script will output a pseudo-random subset of the subject file (I<i.e.>, arg3), such that its distribution matches T's density as closely as possible.

=head1 DESCRIPTION

Given distinct B<"subject" (S)> and a B<"target" (T)> distributions, this script attempts to mimic T's density (I<i.e.>, its shape) by pseudo-randomly sampling from S's population.

Warning: The script attempts to match T's density only, not its population size.

IMPORTANT NOTES

=over

=item B<"Pseudo-randomness">

Items from S's population are randomly selected within bins, not within the entire population, hence the "pseudo" prefix

=item B<Number of bins to choose>

Usually the more, the better.

=item B<Log-transform>

Binning into log10-transformed is highly recommended I<e.g.> for matching FPKM/RPKM distributions (see transform option).

=item B<Re-iteration>

It might be necessary to call the script several times sequentially (I<i.e.> input -> output1; output1 -> output2; output2 -> output3, etc., where "->" denotes a matchDistribution call) until reaching an optimum.
This is what the accompanying B<matchDistributionLoop.sh> script does (see below).

=back

=head1 RE-ITERATIONS

Use B<matchDistributionLoop.sh> and B<matchDistributionKStest.r>. Both scripts need to be in your $PATH.

B<Usage>: C<< matchDistributionLoop.sh <passes> <doKolmogorov-Smirnov> <target> <subject> <bins> <breakIfKSTest> >>

Where:

=over

=item B<passes> (int): Maximum number of passes to perform

=item B<doKolmogorov-Smirnov> (0|1 boolean): Toggle do KS test on resulting distributions after each pass, and print p-value. This will call B<matchDistributionKStest.r> (courtesy of Andres Lanzos, CRG).

=item B<target> (string): Path to file containing T's values.

=item B<subject> (string): Path to tab-separated file containing S's identifiers and values.

=item B<bins> (int): Number of bins to split the distributions into.

=item B<breakIfKSTest> (0|1 boolean): break loop if KS test gives p>0.05 (I<i.e.>, before reaching the maximum number of passes).

=back

=head1 DEPENDENCIES

CPAN: List::Util 'shuffle'

=head1 AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

=cut

my $transform;
my $verbose='';
GetOptions (
        'transform=s' => \$transform,
        'verbose' => \$verbose
        );

our $pseudocount=0.001;

if(defined $transform){
	unless ($transform eq 'log10'){
		$message_text="invalid transform value (must be 'log10')\n";
		pod2usage( { -message => $message_text,
        		     -exitval => $exit_status  ,
               		-verbose => $verbose_level,
               -output  => $filehandle } );
		  }
	print STDERR "Working on $transform - transformed values (with pseudocount = $pseudocount). Output values will be the original, non-transformed ones, though.\n";
}
else{
	$transform='no';
}

# http://stackoverflow.com/questions/1915616/how-can-i-elegantly-call-a-perl-subroutine-whose-name-is-held-in-a-variable
my %transform_subs = (no => \&no,
                      log10 => \&log10
                      );

my $distribToMatch=$ARGV[0]; # one column with values
my $bins=$ARGV[1]; # number of bins to separate the distrib into
my $targetObjectsToSubsampleFrom=$ARGV[2]; # two columns, e.g. column 1 is gene_id, column2 is RPKM
my @origvalues=();

if($#ARGV!=2){
	$message_text="Wrong number of arguments.\n";
pod2usage( { -message => $message_text,
        		     -exitval => $exit_status  ,
               		-verbose => $verbose_level,
               -output  => $filehandle } );
		  }

open D, "$distribToMatch" or die $!;

print STDERR "Reading target distribution to mimic...\n";
while(<D>){
	chomp;
	push(@origvalues, $transform_subs{$transform}->($_));
		if ($.%1000000 == 0){
                print STDERR "\tProcessed $. lines\n";
        }

}

#sort array numerically
my @sorted = sort { $a <=> $b } @origvalues;
@origvalues=@sorted;
@sorted=();
close D;

my $min=$origvalues[0];
my $max=$origvalues[$#origvalues];
my $range=$max-$min;
my $binRange=$range/$bins;

print STDERR "Min: $min, Max: $max. Range: $range. Bin range: $binRange\n";

#calculate what fraction of the total each bin represents

my %binSizes=();
for (my $i=$min; $i<=$max; $i+=$binRange){
	my $count=0;
	foreach my $j (@origvalues){
		if($j>=$i && $j<$i+$binRange){
			$count++
		}
		elsif($j>$i+$binRange){
			last;
		}
	}
	$binSizes{$i}=$count/($#origvalues+1);
}
print STDERR "Done...\n";

my %targetIdsToValues=();
my %targetsList=();
print STDERR "Reading subject dataset to sample from...\n";
open T, "$targetObjectsToSubsampleFrom" or die $!;
my %binnedPop=();
my $countValuesWithinRange=0;
while(<T>){
	chomp;
	my @line=split "\t";
	$targetIdsToValues{$line[0]}=$line[1];
	#populate bins
	foreach my $bin (keys %binSizes){
		if($transform_subs{$transform}->($line[1])>=$bin && $transform_subs{$transform}->($line[1])<$bin+$binRange){
			push(@{$binnedPop{$bin}}, $line[0]);
			$countValuesWithinRange++;
			last;
		}

	}
	if ($.%1000000 == 0){
                print STDERR "\tProcessed $. lines\n";
        }
}
close T;
print STDERR "Done...\n";


print STDERR "Found $countValuesWithinRange values within range in $targetObjectsToSubsampleFrom.\n" if($verbose);

foreach my $bin (keys %binnedPop){
	my $size=$#{$binnedPop{$bin}}+1;
	print STDERR "\nbin $bin: N= $size\n" if($verbose);
	print STDERR " Fraction of input: ".$size/$countValuesWithinRange."\n" if($verbose);
	print STDERR " Desired fraction in output: ".$binSizes{$bin}."\n" if($verbose);
	my $numberOfItemsToPick=int($binSizes{$bin}*$countValuesWithinRange);
	$numberOfItemsToPick=$size if ($size<$numberOfItemsToPick);
	print STDERR " Will try to pick $numberOfItemsToPick items at random.\n" if($verbose);
	print STDERR "###### WARNING bin $bin will be of size 0!! (No values available in input)\n" if ($size==0 && $verbose);
	my @shuffled = shuffle(@{$binnedPop{$bin}});
	my @ids=splice(@shuffled, 0, $numberOfItemsToPick);
	foreach my $id (@ids){
		print "$id\t$targetIdsToValues{$id}\n";
	}
}
print STDERR "Done!\n";

sub log10 {
        my $n = shift;
        #add pseudocount
        return log($n+$pseudocount)/log(10);
    }
sub no {
        my $n = shift;
        return $n;
    }