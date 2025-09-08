use strict;
use warnings;

use Pod::Usage;

my $message_text  = "Error\n";
my $exit_status   = 2;          ## The exit status to use
my $verbose_level = 99;          ## The verbose level to use
my $filehandle    = \*STDERR;   ## The filehandle to write to
my $sections = "NAME|SYNOPSIS|DESCRIPTION";



=head1 NAME

hashToGff

=head1 SYNOPSIS

Converts a multidimensional hash, as produced by gffToHash.pm, into an array containing one GFF line (with endline "\n" included) per element.

B<Usage>: C<< @outGff = hashToGff(<gff_hashref>) >>

=cut


sub hashToGff{
	my %gffHash=%{$_[0]};
	my @outGff=();
	foreach my $id (keys %gffHash){
  		for (my $j=0; $j<=$#{$gffHash{$id}};$j++){
			my @outGffLine=();
			my @attrs;
			for (my $i=0; $i<=7; $i++){ #processing the first 8 GFF fields
			  push(@outGffLine, $gffHash{$id}[$j][$i]);
			}
			my @sortedKeys=sort keys %{${$gffHash{$id}[$j]}[8]};
			my @resortedKeys=();
			my @mainKeys=(); # gene_id  and transcript_id, to be prepended in that order to @resortedKeys
			for (my $i=0; $i<=$#sortedKeys; $i++){
				if($sortedKeys[$i] eq 'transcript_id'){
					push(@mainKeys, $sortedKeys[$i]);
				}
				elsif($sortedKeys[$i] eq 'gene_id'){
					unshift(@mainKeys, $sortedKeys[$i])
				}
				else{
					push(@resortedKeys, $sortedKeys[$i]);
				}
			}
			unshift(@resortedKeys, @mainKeys);

			foreach my $key (@resortedKeys){ #processing GFF attributes (9th field)
			    push(@attrs, $key." \"${${$gffHash{$id}[$j]}[8]}{$key}\";");
			}

			my $attr=join(" ", @attrs);
			push(@outGffLine, $attr);
			push(@outGff, join("\t", @outGffLine)."\n");
		}
	}
	return @outGff;
}



1;
