
if config["LORDEC_CORRECT"]:
	rule splitLrFastqs:
		input: FQPATH + "{techname}_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQPATH + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
		output: temp(FQPATH + "splitFasta/{techname}_{capDesign}_{sizeFrac}_{split}.fasta" if config["DEMULTIPLEX"] else temp(FQPATH + "splitFasta/{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.fasta"))
		shell:
			'''
let splitInto={splitFastqsInto}+1
uuid=$(uuidgen)
zcat {input} | fastq2fasta.pl|FastaToTbl > {config[TMPDIR]}/$uuid
split -a 4 --additional-suffix .tbl -d -n l/$splitInto {config[TMPDIR]}/$uuid {config[TMPDIR]}/${{uuid}}_
cat {config[TMPDIR]}/${{uuid}}_{wildcards.split}.tbl | TblToFasta > {output}
			'''

if config["LORDEC_CORRECT"]:
	rule lordecCorrectLr:
		input:
			lr=FQPATH + "splitFasta/{techname}_{capDesign}_{sizeFrac}_{split}.fasta"  if config["DEMULTIPLEX"] else FQPATH + "splitFasta/{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.fasta",
			sr1=config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_1.fastq.gz" if config["DEMULTIPLEX"] else lambda wildcards: expand(config["MATCHED_HISEQ_PATH"] + "{techname}_{capDesign}.{barcodes}_1.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes),
			sr2=config["MATCHED_HISEQ_PATH"] + "hiSeq_{capDesign}_2.fastq.gz" if config["DEMULTIPLEX"] else lambda wildcards: expand(config["MATCHED_HISEQ_PATH"] + "{techname}_{capDesign}.{barcodes}_2.fastq.gz", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, barcodes=wildcards.barcodes)
		output: temp(FQPATH + "splitFasta/corr/{techname}_{capDesign}_{sizeFrac}_{split}.Corr" +  lastK + ".fasta" if config["DEMULTIPLEX"] else FQPATH + "splitFasta/corr/{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.Corr" +  lastK + ".fasta")
		threads: 6
		shell:
			'''
kvalues={graph_kmers_string}
for i in "${{!kvalues[@]}}"; do
kval=${{kvalues[$i]}}
echoerr "kmer $kval"
echoerr "Building SR graph"
lordec-build-SR-graph -T {threads} -O {config[TMPDIR]}/ -2 {input.sr1},{input.sr2} -k $kval -s {solid_kmer_abundance_threshold} -g {config[TMPDIR]}/hiSeq_{wildcards.capDesign}_k${{kval}}_s{solid_kmer_abundance_threshold}.h5

outFileBn=$(basename {output} .Corr{lastK}.fasta).Corr${{kval}}.fasta
outFileDn=$(dirname {output})
outFile=$outFileDn"/"$outFileBn

if [ $i -gt 0 ]; then
previousKval=${{kvalues[$i-1]}}

previousFileBn=$(basename {output} .Corr{lastK}.fasta).Corr${{previousKval}}.fasta
fileDn=$(dirname {output})
previousFile=$fileDn"/"$previousFileBn


lordec-correct -T {threads} -O {config[TMPDIR]}/ -2 {config[TMPDIR]}/hiSeq_{wildcards.capDesign} -k $kval -s {solid_kmer_abundance_threshold} -i $previousFile -o $outFile


else
lordec-correct -T {threads} -O {config[TMPDIR]}/ -2 {config[TMPDIR]}/hiSeq_{wildcards.capDesign} -k $kval -s {solid_kmer_abundance_threshold} -i {input.lr} -o $outFile
fi
rm -f {config[TMPDIR]}/hiSeq_{wildcards.capDesign}_k${{kval}}_s{solid_kmer_abundance_threshold}.h5
done

			'''

if config["LORDEC_CORRECT"]:
	rule mergeCorrLrFastas:
		input: lambda wildcards: expand(FQPATH + "splitFasta/corr/{techname}_{capDesign}_{sizeFrac}_{split}.Corr" + lastK + ".fasta", techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, split=splitFasta) if config["DEMULTIPLEX"] else expand(FQPATH + "splitFasta/corr/{techname}_{capDesign}_{sizeFrac}.{barcodes}_{split}.Corr" + lastK + ".fasta", filtered_product, techname=wildcards.techname, capDesign=wildcards.capDesign, sizeFrac=wildcards.sizeFrac, barcodes=wildcards.barcodes, split=splitFasta)
		output: FQ_CORR_PATH + "{techname}Corr" + lastK + "_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr" + lastK + "_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
		shell:
			'''
cat {input} | FastaToTsv | tsv2fastq.pl |gzip > {output}

			'''

rule linkOriginalFastqs:
	input: FQPATH + "{techname}_{capDesign}_{sizeFrac}.fastq.gz"  if config["DEMULTIPLEX"] else FQPATH + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "{techname}Corr" + FINALCORRECTIONLEVELS[0] + "_{capDesign}_{sizeFrac}.fastq.gz" if config["DEMULTIPLEX"] else FQ_CORR_PATH + "{techname}Corr" + FINALCORRECTIONLEVELS[0] + "_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	shell:
		'''
ln -sr {input} {output}
		'''
