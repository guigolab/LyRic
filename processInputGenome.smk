


rule sortIndexGenome:
	input: config["GENOMESDIR"] +"{genome}.fa"
	output: 
		sorted=config["GENOMESDIR"] +"{genome}.sorted.fa",
		bioperlindex=config["GENOMESDIR"] +"{genome}.sorted.fa.index"
	conda: "envs/perl_env.yml"
	shell:
		'''
uuid=$(uuidgen)
#check for duplicate sequences:
count=$(cat {input} | fgrep ">" | sort|uniq -d | wc -l)
if [ $count -gt 0 ]; then echo "$count duplicate sequence IDs found"; exit 1; fi

FastaToTbl {input} | sort -T {TMPDIR} -k1,1 | TblToFasta > {TMPDIR}/$uuid 
mv {TMPDIR}/$uuid {output.sorted}
perl -e 'use Bio::DB::Fasta; my $chrdb = Bio::DB::Fasta->new("{output.sorted}");'

		'''

rule makeGenomeFile:
	input: config["GENOMESDIR"] +"{genome}.sorted.fa"
	output: config["GENOMESDIR"] +"{genome}.sorted.genome"
	shell:
		'''
 uuid=$(uuidgen)
FastaToTbl {input} | awk '{{print $1"\\t"length($2)}}' | sort -k1,1 > {TMPDIR}/$uuid
mv {TMPDIR}/$uuid {output}
		'''

rule makeGencodePartition:
	input:
		gtf=lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]],
		genome = lambda wildcards: config["GENOMESDIR"] + CAPDESIGNTOGENOME[wildcards.capDesign] + ".sorted.genome"
	conda: "envs/xtools_env.yml"
	output: "output/annotations/{capDesign}.partition.gff"
	shell:
		'''

partitionAnnotation.sh {input.gtf} {input.genome} | sortgff > {output}

# QC:
genomeSize=$(cat {input.genome} | cut -f2|sum.sh)
testSize=$(cat {output} | awk '{{print ($5-$4)+1}}' | sum.sh)
if [ $testSize -ne $genomeSize ]; then
echoerr "ERROR: sum of feature sizes in output gff is not equal to genome size. The output is probably bogus."
exit 1
fi
		'''


rule makeSirvGff:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/{capDesign}.SIRVs.gff"
	shell:
		'''
cat {input} | awk '$1 ~ /SIRV/' | sortgff > {output}

		'''

rule makeGencodeBed:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/{capDesign}.bed"
	shell:
		'''
cat {input} | gff2bed_full.pl - |sortbed > {output}
		'''

rule simplifyGencode:
	input: lambda wildcards: GENOMETOANNOTGTF[CAPDESIGNTOGENOME[wildcards.capDesign]]
	output: "output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	shell:
		'''
uuidTmpOut=$(uuidgen)
cat {input}  | simplifyGencodeGeneTypes.pl - | sort -T {TMPDIR}	 -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}
		'''

rule collapseGencode:
	input: "output/annotations/simplified/{capDesign}.gencode.simplified_biotypes.gtf"
	output: "output/annotations/simplified/{capDesign}.gencode.collapsed.simplified_biotypes.gtf"
	threads:1
	conda: "envs/xtools_env.yml"
	shell:
		'''
uuid=$(uuidgen)
uuidTmpOut=$(uuidgen)
cat {input} | skipcomments | sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  | tmerge --exonOverhangTolerance {ExonOverhangTolerance} - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuid
uuidL=$(uuidgen)

bedtools intersect -s -wao -a {TMPDIR}/$uuid -b {TMPDIR}/$uuid | buildLoci.pl - |sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n	> {TMPDIR}/$uuidL
mergeToRef.pl {input} {TMPDIR}/$uuidL | sort -T {TMPDIR}  -k1,1 -k4,4n -k5,5n  > {TMPDIR}/$uuidTmpOut
mv {TMPDIR}/$uuidTmpOut {output}

		'''


