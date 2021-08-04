
rule linkOriginalFastqs:
	input: LR_FASTQDIR + "{techname}_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	output: FQ_CORR_PATH + "{techname}Corr" + FINALCORRECTIONLEVELS[0] + "_{capDesign}_{sizeFrac}_{barcodes}.fastq.gz"
	shell:
		'''
ln -sr {input} {output}
		'''
