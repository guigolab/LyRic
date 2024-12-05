#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

no warnings 'recursion';

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";

=head1 NAME

buildLoci

=head1 SYNOPSIS

A utility to build gene loci (i.e., sets of overlapping transcripts) out of a transcript set.

B<Usage example>:

C<< bedtools intersect -s -wao -a inGTF -b inGTF | buildLoci.pl - > test.loci.gff >>

=head2 INPUT

The input file is provided as first argument to the script. It consists in two ("left" and "right") GTF records per line, separated by a C<< tab >>. This is typically the standard output produced by C<< bedtools intersect -wao -a inGTF -b inGTF >> (C<< inGTF >> being a single GTF file).

The flexibility of C<< bedtools >> allows the user to build gene loci based on whatever definition they need, e.g., with or without respect to genomic strand (see C<< bedtools >>'s C<< -s >> option).

=head2 OPTIONS

=over

=item B<keepGeneid> = If set, any C<< gene_id >> values present in the input will be kept in the output, under attribute C<< gene_id_bkp >>.


=item B<locPrefix> (string) = When set, this parameter's value will be prepended to all gene_id values in the output (in the form of C<< <locPrefix>LOC_XXXXXXXXXX >>)

=back

=head2 OUTPUT

One GTF line per unique "left" GTF record in the input, with a supplementary C<< gene_id >> attribute (in the form of C<< LOC_XXXXXXXXXX >>) appended to its 9th field.

Any gene_id value present in the input will be overwritten, except if the C<< --keepGeneid >> option is used.

=head1 DESCRIPTION

Any pair of GTF records present within a line of input is assumed to represent overlapping features. Gene loci are then built based on these overlaps, i.e. the <transcript_id>s of both records are assigned the same arbitrary <gene_id> value in the output.

=head1 DEPENDENCIES

Although not strictly necessary, L<BEDTools|https://github.com/arq5x/bedtools2> is recommended, solely to provide input to C<< buildLoci.pl >>.

=head1 AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

=cut

my $keep_gene_id_backup='';
my $locPrefix='';
GetOptions ('keepGeneid' => \$keep_gene_id_backup,
            'locPrefix=s' => \$locPrefix
            )
or pod2usage( { -message => "Error in command line arguments",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );

unless (defined $ARGV[0]){
	pod2usage( { -message => "Error in command line arguments: no input provided.",
        		  -exitval => $exit_status  ,
            		-verbose => $verbose_level,
               -output  => $filehandle } );
}

my %gffOut;
my %transcript_to_transcripts=();
open GFFINT, "$ARGV[0]" or die $!;
print STDERR "Parsing input...\n";
while(<GFFINT>){
	chomp;
	my $line=$_;
	my @line=split("\t", $line);
	if($keep_gene_id_backup){
		$line[8]=~s/gene_id (\"\S+\")/gene_id_bkp $1/g;
	}
	else{
		$line[8]=~s/gene_id \"\S+\";//g;

	}
	my @gffRecord=@line[0..8];
	$gffOut{join("\t",@gffRecord)}=1;
	if($line=~/transcript_id \"(\S+)\";.+transcript_id \"(\S+)\";/){
		my $tr1=$1;
		my $tr2=$2;
		$transcript_to_transcripts{$tr1}{$tr2}=1;
	}
	else{
		die "Line must contain two transcript_id attributes. Offending line: $. . Can't continue.\n"
	}
}
close GFFINT;
print STDERR "Parsing input DONE.\n";
my %transcript_id_to_locus_id=();
my $locusNumber=0;
foreach my $tr1 (keys %transcript_to_transcripts){
	searchOverlap($tr1);
	$locusNumber++;
}
foreach my $gffRecord ( keys %gffOut){
	my @gffRecord=split("\t", $gffRecord);
	$gffRecord[8]=~/transcript_id \"(\S+)\"/;
	my $tr=$1;
	$gffRecord[8].=" gene_id \"".makeLocusId($transcript_id_to_locus_id{$tr})."\";";
	print join ("\t", @gffRecord)."\n";
}
sub searchOverlap{
	my $tr1=$_[0];
	unless(exists $transcript_id_to_locus_id{$tr1}){
		$transcript_id_to_locus_id{$tr1}=$locusNumber;
	}
	foreach my $tr2 (keys %{$transcript_to_transcripts{$tr1}})
	{
		unless( exists $transcript_id_to_locus_id{$tr2}){
			searchOverlap($tr2);
		}
	}
}

sub makeLocusId{
	my $id=$_[0];
	my @newId=split("", $id);
	#my @prefix=split("", $locPrefix);
	my @prepend=(join("",split("", $locPrefix)),'L','O','C','_');
	my $totalLength=($#prepend+13)-length($id);
	#print STDERR "$#prepend $tmp\n";
	for (my $i=$#prepend+1;$i<$totalLength; $i++){
	#	print STDERR "\t$i\n";
		$prepend[$i]=0;
	}
	unshift(@newId, @prepend);
	return(join("",@newId))
}
