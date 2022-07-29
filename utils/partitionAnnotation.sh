#!/bin/bash

gtf=$1
genome=$2
uuid=$(uuidgen)

cat $genome | cut -f1 | sort |uniq > $uuid.chr.list

zcat -f $gtf | fgrep -w -f $uuid.chr.list | sortgff > $uuid.ingenomegtf

# projected exons:
zcat -f $uuid.ingenomegtf |  awk '$3=="exon"' | sortgff | bedtools merge -i stdin > $uuid.exons.nostrand.bed
# projected CDSs
zcat -f $uuid.ingenomegtf |  awk '$3=="CDS"' | fgrep "transcript_type \"protein_coding\";" | sortgff | bedtools merge -i stdin > $uuid.CDS.nostrand.bed

# projected exons of coding transcripts:
zcat -f $uuid.ingenomegtf |  awk '$3=="exon"' | fgrep "transcript_type \"protein_coding\";"| sortgff | bedtools merge -i stdin > $uuid.exonsOfCodingTranscripts.nostrand.bed

# projected exons of pseudogenes
zcat -f $uuid.ingenomegtf |  awk '$3=="exon"' | grep -P 'transcript_type "\S+seudogene\S*";'| sortgff | bedtools merge -i stdin | bedtools subtract -a stdin -b $uuid.exonsOfCodingTranscripts.nostrand.bed > $uuid.exonsOfPseudogenes.nostrand.bed


# projected UTRs of coding transcripts:
bedtools subtract -a $uuid.exonsOfCodingTranscripts.nostrand.bed -b $uuid.CDS.nostrand.bed | sortbed | bedtools merge -i stdin > $uuid.UTRs.nostrand.bed


#projected exons of noncoding transcripts:
bedtools subtract -a $uuid.exons.nostrand.bed -b $uuid.exonsOfCodingTranscripts.nostrand.bed | sortbed | bedtools merge -i stdin | bedtools subtract -a stdin -b $uuid.exonsOfPseudogenes.nostrand.bed > $uuid.exonsOfNonCodingTranscripts.nostrand.bed

# introns that are never exons
zcat -f $uuid.ingenomegtf |  awk '$3=="exon"' | makeIntrons.pl - |sortgff| bedtools merge -i stdin | bedtools subtract -a stdin -b $uuid.exons.nostrand.bed > $uuid.introns.nostrand.bed

# genic regions
cat $uuid.exons.nostrand.bed $uuid.introns.nostrand.bed | sortbed | bedtools merge -i stdin > $uuid.genic.nostrand.bed

#intergenic regions

bedtools complement -i $uuid.genic.nostrand.bed -g $genome | sortbed | bedtools merge -i stdin > $uuid.intergenic.nostrand.bed

#make GFF

cat $uuid.intergenic.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"intergenic\";"}' > $uuid.out.gff
cat $uuid.introns.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"intron\";"}' >> $uuid.out.gff
cat $uuid.CDS.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"CDS\";"}' >> $uuid.out.gff
cat $uuid.UTRs.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"UTR\";"}' >> $uuid.out.gff
cat $uuid.exonsOfNonCodingTranscripts.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"exonOfNCT\";"}' >> $uuid.out.gff
cat $uuid.exonsOfPseudogenes.nostrand.bed | awk '{print $0"\t"$1"_"$2+1"_"$3"\t0\t."}' | bed2gff.pl | awk '{print $0" region_flag \"exonOfPseudo\";"}' >> $uuid.out.gff

cat $uuid.out.gff | sortgff 

###################
##### QC ##########
###################
genomeSize=$(cat $genome | cut -f2|sum.sh)
intergenicSize=$(cat $uuid.intergenic.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
genicSize=$(cat $uuid.genic.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
test=$((intergenicSize+genicSize))
if [ $test -ne $genomeSize ]; then
echoerr "ERROR: sum of (genic + intergenic) sizes is not equal to genome size. The output is probably bogus."
rm -f $uuid.*
exit 1
fi
intronsSize=$(cat $uuid.introns.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
cdsSize=$(cat $uuid.CDS.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
utrsSize=$(cat $uuid.UTRs.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
enctSize=$(cat $uuid.exonsOfNonCodingTranscripts.nostrand.bed|awk '{ print $3-$2}'|sum.sh)
epsSize=$(cat $uuid.exonsOfPseudogenes.nostrand.bed |awk '{ print $3-$2}'|sum.sh)
test2=$((intronsSize+cdsSize+utrsSize+enctSize+epsSize))
if [ $test2 -ne $genicSize ]; then
echoerr "ERROR: sum of (introns + cds + utrs + exonsOfNonCodingTranscripts + exonsOfPseudogenes) sizes is not equal to genic size. The output is probably bogus."
rm -f $uuid.*
exit 1
fi
rm -f $uuid.*
