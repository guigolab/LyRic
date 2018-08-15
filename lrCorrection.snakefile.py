
if config["LORDEC_CORRECT"]:
	rule splitLrFastqs:
		input: FQPATH + "{techname}_{capDesign}_{sizeFrac}.fastq.gz"
		output: temp(FQPATH + "splitFasta/" + "{techname}_{capDesign}_{sizeFrac}_{split}.fasta" if config["DEMULTIPLEX"] else "splitFasta/" + "{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.fasta")
		shell:
			'''
let splitInto={splitFastqsInto}+1
zcat {input} | fastq2fasta.pl|FastaToTbl > $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.tbl
split -a 4 --additional-suffix .tbl -d -n l/$splitInto $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}.tbl  $TMPDIR/{wildcards.techname}_{wildcards.capDesign}_{wildcards.sizeFrac}_
cat $TMPDIR/$(basename {output} .fasta).tbl | TblToFasta > {output}
			'''

if config["LORDEC_CORRECT"]:
	rule lordecCorrectLr:
		input:
			lr=FQPATH + "splitFasta/" + "{techname}_{capDesign}_{sizeFrac}_{split}.fasta"  if config["DEMULTIPLEX"] else "splitFasta/" + "{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.fasta",
			sr1="fastqs/" + "hiSeq_{capDesign}_1.fastq.gz",
			sr2="fastqs/" + "hiSeq_{capDesign}_2.fastq.gz"
		output: temp(FQPATH + "splitFasta/corr/" + "{techname}_{capDesign}_{sizeFrac}_{split}.Corr" +  lastK + ".fasta" if config["DEMULTIPLEX"] else FQPATH + "splitFasta/corr/" + "{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.Corr" +  lastK + ".fasta")
		threads: 6
		shell:
			'''
kvalues={graph_kmers_string}
for i in "${{!kvalues[@]}}"; do
kval=${{kvalues[$i]}}
echoerr "kmer $kval"
echoerr "Building SR graph"
lordec-build-SR-graph -T {threads} -O $TMPDIR/ -2 {input.sr1},{input.sr2} -k $kval -s {solid_kmer_abundance_threshold} -g $TMPDIR/hiSeq_{wildcards.capDesign}_k${{kval}}_s{solid_kmer_abundance_threshold}.h5

outFileBn=$(basename {output} .Corr{lastK}.fasta).Corr${{kval}}.fasta
outFileDn=$(dirname {output})
outFile=$outFileDn"/"$outFileBn

if [ $i -gt 0 ]; then
previousKval=${{kvalues[$i-1]}}

previousFileBn=$(basename {output} .Corr{lastK}.fasta).Corr${{previousKval}}.fasta
fileDn=$(dirname {output})
previousFile=$fileDn"/"$previousFileBn


lordec-correct -T {threads} -O $TMPDIR/ -2 $TMPDIR/hiSeq_{wildcards.capDesign} -k $kval -s {solid_kmer_abundance_threshold} -i $previousFile -o $outFile


else
lordec-correct -T {threads} -O $TMPDIR/ -2 $TMPDIR/hiSeq_{wildcards.capDesign} -k $kval -s {solid_kmer_abundance_threshold} -i {input.lr} -o $outFile
fi
rm -f $TMPDIR/hiSeq_{wildcards.capDesign}_k${{kval}}_s{solid_kmer_abundance_threshold}.h5
done

			'''

if config["LORDEC_CORRECT"]:
	rule mergeCorrLrFastas:
		input: lambda wildcards: expand(FQPATH + "splitFasta/corr/" + "{techname}_{capDesign}_{sizeFrac}_{split}.Corr" + lastK + ".fasta", techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, split=splitFasta) if config["DEMULTIPLEX"] else expand(FQPATH + "splitFasta/corr/" + "{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.Corr" + lastK + ".fasta", techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=BARCODES, split=splitFasta)
		output: FQ_CORR_PATH + "{techname}Corr" + lastK + "_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr" + lastK + "_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
		shell:
			'''
cat {input} | FastaToTsv | tsv2fastq.pl |gzip > {output}

			'''

rule linkOriginalFastqs:
	input: FQPATH + "{techname}_{capDesign}_{sizeFrac}.fastq.gz"  if config["DEMULTIPLEX"] else FQPATH + "{techname}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "{techname}Corr" + FINALCORRECTIONLEVELS[0] + "_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr{corrLevel}_{capDesign}_{sizeFrac}.{barcodes}.fastq.gz"
	shell:
		'''
ln -sr {input} {output}
		'''
