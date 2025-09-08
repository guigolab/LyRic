#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use FindBin;    # find present script
use lib "$FindBin::Bin";
use gffToHash;
use hashToGff;


my @biotypePriority=("SIRV","ERCC", "pseudogene", "protein_coding", "rRNA", "ribozyme", "lncRNA", "tRNA", "vaultRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA","ribozyme","scaRNA","scRNA","snoRNA","snRNA","sRNA","TEC","artifact","vault_RNA");


my $refGtf=$ARGV[0]; #refenrece gencode gtf, e.g. (simplified for readability):
#chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "pseudogene";

my $mergedLoci=$ARGV[1];
#chr17   tmerge  exon    35147808        35147910        .       +       .        transcript_id "TM_000000110823"; contains "cls.TM_000000036661,cls.TM_000000036663,ENST00000456328.2"; gene_id "LOC_000000020312";

open REF, "$refGtf" or die $!;

my %ref_transcript_id_to_gene_id=();
my %ref_gene_id_to_gene_type=();

while(<REF>){
	my $transcript_id;
	my $gene_id;
	my $gene_type;
	if($_=~/transcript_id \"(\S+)\";/){
		$transcript_id=$1;
	}
	if($_=~/gene_id \"(\S+)\";/){
		$gene_id=$1;
	}
	if($_=~/gene_type \"(\S+)\";/){
		$gene_type=$1;
	}
	if(defined $transcript_id){
		if(defined $gene_id){
			$ref_transcript_id_to_gene_id{$transcript_id}=$gene_id;
			if(defined $gene_type){
				$ref_gene_id_to_gene_type{$gene_id}=$gene_type;
			}
			else{
				die "No gene_type attribute found for $gene_id, can't continue.\n";
			}
		}
	}

}

close REF;

my %gffHash=gffToHash($mergedLoci, 'transcript_id', 0,);

#print Dumper \%gffHash;

my %merged_gene_idToRefgene_id=();
#my %merged_gene_idToRefgene_type=();


#map merged loci IDs to reference gene_id's if they exist
foreach my $tx (keys %gffHash){
	foreach my $exon (@{$gffHash{$tx}}){
		#print ${$exon}[8]{'contains'}."\n";
		my $contains= ${$exon}[8]{'contains'};
		my $gene_id=${$exon}[8]{'gene_id'};
		#print "$tx,$contains\n";
		my @cont=split(",", $contains);
		foreach my $contTx (@cont){
			if(exists $ref_transcript_id_to_gene_id{$contTx}){
				my $refgene_id=$ref_transcript_id_to_gene_id{$contTx};
				$merged_gene_idToRefgene_id{$gene_id}{$refgene_id}=1;
			}

		}
	}

}

foreach my $tx (keys %gffHash){
	foreach my $exon (@{$gffHash{$tx}}){
		my $gene_type;
		my $gene_ref_status;
		my $gene_id=${$exon}[8]{'gene_id'};
		if(exists ($merged_gene_idToRefgene_id{$gene_id})){
			my @ref_gene_ids=();
			foreach my $ref_gene_id (sort keys %{$merged_gene_idToRefgene_id{$gene_id}}){
				push(@ref_gene_ids, $ref_gene_id)
			}
			${$exon}[8]{'ref_gene_ids'}=join(",", @ref_gene_ids);

			my $found=0;
			foreach my $registered_genetype (@biotypePriority){
				if($found == 0){
					foreach my $refGeneId (@ref_gene_ids){
						if($registered_genetype eq $ref_gene_id_to_gene_type{$refGeneId}){
							${$exon}[8]{'gene_id'}=$refGeneId;
							${$exon}[8]{'gene_type'}=$ref_gene_id_to_gene_type{$refGeneId};
							$found=1;
							last;
						}
					}
				}
			}
			if($found==1){
				$gene_ref_status="known";
			}
			else{
				die "Unkown gene_type for $gene_id. ".'Edit \@biotypePriority ?\n';
			}

		}
		else{
			$gene_ref_status="novel";
			${$exon}[8]{'gene_type'}="lncRNA";
		}
		${$exon}[8]{'gene_ref_status'}=$gene_ref_status;


	}
}


my @outGff=hashToGff(\%gffHash);
print join("", @outGff);
