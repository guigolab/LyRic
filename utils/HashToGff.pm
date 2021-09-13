use strict;
use warnings;

#takes a GFF as input and returns a multidimensional hash looking like this (as reported by Data::Dumper):
# $VAR1 = {
#         'chr5_78623031_78623031_+' => [
#                                          [
#                                            'chr5',
#                                            'RnaSeq',
#                                            'polyAsite',
#                                            '78623031',
#                                            '78623031',
#                                            '.',
#                                            '+',
#                                            '.',
#                                            {
#                                              'ln' => '80',
#                                              'id' => 'chr5_78623031_78623031_+'
#                                            }
#                                          ],
#                                          [
#                                           'chr5',
#                                           'RnaSeq',
#                                           'polyAsite',
#                                           '78623031',
#                                           '78623031',
#                                           '.',
#                                           '+',
#                                           '.',
#                                           {
#                                             'ln' => '81',
#                                             'id' => 'chr5_78623031_78623031_+'
#                                           }
#                                          ]
#                                          ],


# optionally (if $_[2] ==1), the whole (chomp'ed) gff line can be kept in the last array element
# optionally (if $_[3] == whatever), only 'whatever' gff records (i.e. 'whatever' in 3rd column of gff) will be output, e.g. 'exon'


sub gffToHash{
	die "Wrong nr of args: @_ .\n" if ($#_<1 || $_[0] eq '' || $_[1] eq '');
	my $keepWholeRecordInLastArrayElement=0;
	if (exists $_[2]){
		$keepWholeRecordInLastArrayElement=$_[2];
	}
	my $onlyGffRecords=undef;
	if (exists $_[3]){
		$onlyGffRecords=$_[3];
	}
	open GFF, "$_[0]" or die $!;
	my %gff=();
#	my %existsElementId=(); #to check if ElementId is really unique within the file
	my $identifyingElementId=$_[1]; # key used to uniquely identify gff records (e.g. exon_id or transcript_id attribute).
	while (<GFF>){
		next if ($_=~/^#/);

		chomp;
		my $untouchedLine=$_;
		my $line=$untouchedLine;
		$line=~s/\"//g;
		$line=~s/;//g;
		my @line=split("\t", $line);
		die "Malformed line at line $. . This doesn't look like a proper 9-field GFF record, dying.\n" unless ($#line==8);
		if(defined $onlyGffRecords){
			next unless ($line[2] eq $onlyGffRecords)
		}
		my $elementId=undef;
		if($line=~/$identifyingElementId (\S+)/){
			$elementId=$1;
#			if(exists $existsElementId{$elementId}){
#				warn "$identifyingElementId '$elementId' is not unique within the file. Only its last occurrence will be returned.\n";
#			}
#			$existsElementId{$elementId}=1;
		}
		else{
			warn "Line $. skipped (contains no '$identifyingElementId' attribute.\n";
			next;
		}
		my @attrs=split(" ", pop(@line));
		die "Odd number of subfields in 9th field of line $.. Dying.\n" unless ($#attrs%2==1);
		push (@{$gff{$elementId}}, \@line);
		#@{$gff{$elementId}}=@line;
		for (my $i=0; $i<=$#attrs;$i=$i+2){
			${${$gff{$elementId}}[$#{$gff{$elementId}}]}[8]{$attrs[$i]}=$attrs[$i+1];
		}
		if($keepWholeRecordInLastArrayElement ==1){
			push (@{$gff{$elementId}[$#{$gff{$elementId}}]}, $untouchedLine);
		}
		#print Dumper \%gff;
	}
	close GFF;
	return %gff;
}

1;
