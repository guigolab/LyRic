#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
#use lib "/users/rg/jlagarde/julien_utils/";


#############################################################################################
# -spliced_coordinates  34 31      30 27    26  22  21  17    16  12   11  7   6    1
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
# +spliced_coordinates  1  4       5  8     9   13  14  18    19  23   24  28  29   34
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
#                       agtc-------cgat-----atcga---tcccg-----gatct----ataca---tgatct
#                       |  |       |  |     |   |   |   |     |   |    |   |   |    |
# chr_coordinates       100|       111|     120 124 128 132   138 142  147 151 155  160
#                          103        114

#                                   gat-----atcga---tcccg-----gat
# +3pprimer_spliced_coords          6                           21   +
# +3pprimer_chr_coords              112                         140  +
#                                   cta-----tagct---agggc-----cta
# +5pprimer_spliced_coords          6                           21   -
# +5pprimer_chr_coords              112                         140  -

#                                   gat-----atcga---tcccg-----gat
# -3pprimer_spliced_coords          29                          14   +
# -3pprimer_chr_coords              112                         140  -
# -5pprimer_spliced_coords          29                          14   -
# -5pprimer_chr_coords              112                         140  +
#############################################################################################




open BEDTSV, "$ARGV[0]" or die $!;
my $previousTranscript_id='';
my $newTranscriptRecordBool=0;
while(<BEDTSV>){
	chomp;
	my @line=split "\t";
	die "Wrong format, line $. of $ARGV[0]. Should be BED12+2 (i.e., the output of 'bedtools coverage -d' on a BED12 file.\n" unless ($#line==13);
	if($line[3] ne $previousTranscript_id){ #new transcript record
		print "New\n"
	}
	else{

	}

}