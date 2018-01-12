import glob
import os.path
from collections import defaultdict
from itertools import product
# # path to dropbox folder where to sync output R plots:
DROPBOX_PLOTS=config["DROPBOX_PLOTSDIR"]
GGPLOT_PUB_QUALITY=config["GGPLOT_PUB_QUALITY"]
#TECHNAME=config["TECHNAME"]
CAPDESIGNTOGENOME=config["capDesignToGenome"]
#mappingDir=config["PB_MAPPINGS"]
sizeFrac_Rpalette=config["SIZEFRACTION_RPALETTE"]
#print(CAPDESIGNTOGENOME)

wildcard_constraints:
 	capDesign = "[^_/]+",
 	sizeFrac = "[^_/]+",
 	techname = "[^_/]+"

# get CAPDESIGNS (capture designs, i.e. Hv1, Hv2, Mv1) and SIZEFRACS (size fractions) variables from FASTQ file names (warning: this will generate duplicate entries):
(TECHNAMES, CAPDESIGNS, SIZEFRACS) = glob_wildcards(config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq")
CAPDESIGNS=set(CAPDESIGNS)
SIZEFRACS=set(SIZEFRACS)
TECHNAMES=set(TECHNAMES)

print (CAPDESIGNS, SIZEFRACS, TECHNAMES)


########################
### make list of forbidden combinations (missing input files)
### inspired by https://stackoverflow.com/questions/41185567/how-to-use-expand-in-snakemake-when-some-particular-combinations-of-wildcards-ar

FORBIDDENCOMBINATIONS = []
AUTHORIZEDCOMBINATIONS = []

for comb in product(TECHNAMES,CAPDESIGNS,SIZEFRACS):
	print (comb)
	if(not os.path.isfile(config["FQPATH"] + comb[0] + "_" + comb[1] + "_" + comb[2] + ".fastq")):
		tup=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]))
		FORBIDDENCOMBINATIONS.append(tup)
	else:
		#tup=frozenset({("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2])})
		tup=(("techname", comb[0]),("capDesign", comb[1]),("sizeFrac", comb[2]))
		AUTHORIZEDCOMBINATIONS.append(tup)
print ("FORBIDDENCOMBINATIONS:")
print (FORBIDDENCOMBINATIONS)
print ("AUTHORIZEDCOMBINATIONS:")
print (AUTHORIZEDCOMBINATIONS)


def filter_combinator(whitelist):
	def filtered_combinator(*args, **kwargs):
		print("args:")
		for a in args:
			print(a)
		print("kwargs:")
		for a in kwargs:
			print(a, kwargs[a])

		for wc_comb in product(*args, **kwargs):
			print ("wc_comb:")
			print (wc_comb)
#			print ("wc_comb 0-3:")
#			print (wc_comb[0:3])

            # Use frozenset instead of tuple
            # in order to accomodate
            # unpredictable wildcard order
            #use only first 2/3 elements of args, which should be techname, capDesign, sizeFrac. We don't care about the rest, as they're not present in input files
			found=False
			for ac in AUTHORIZEDCOMBINATIONS:
#				print("current AUTH:")
#				print (ac)
				if(wc_comb[0:3] == ac):
					print ("SUCCESS")
					found=True
					yield(wc_comb)
					break
			if not found and wc_comb[2][0] != 'sizeFrac':
#if this point is reached it means wc_comb[0:3] was not found. Try with wc_comb[0:2]
				for ac in AUTHORIZEDCOMBINATIONS:
#				print("current AUTH:")
#				print (ac)
					if(wc_comb[0:2] == ac[0:2]):
						print ("SUCCESS")
						found=True
						yield(wc_comb)
						break
#				for currAc in ac:
#					print(currAc)

			# if (frozenset(wc_comb[0:3]) in whitelist):
			# 	print(" OK")
			# 	yield wc_comb
			# # elif (frozenset(wc_comb[0:2]) not in blacklist):
			# # 	print(" OK")
			# # 	yield wc_comb
			# else:
			# 	print(" FORBIDDEN")
	return filtered_combinator

#filtered_product = filter_combinator(FORBIDDENCOMBINATIONS)
filtered_product = filter_combinator(AUTHORIZEDCOMBINATIONS)


rule all:
	input:
#		expand("{techname}_{capDesign}_{sizeFrac}.readlength.tsv", filtered_product,  techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS),
		#expand(expand(config["PLOTSDIR"] + "{techname}_{capDesign}_all.readlength.{{ext}}", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS), ext=config["PLOTFORMATS"])
		expand(config["STATSDATADIR"] + "{techname}_{capDesign}_all.readlength.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
#        expand(expand("results/{{domain}}/{supergroup}/{supergroup_other}/OGSmapping.txt.list.{{cutoff}}.statistics", filtered_product, supergroup=SUPERGROUPS, supergroup_other = SUPERGROUPS_INVERSED), cutoff=CUTOFFS, domain = DOMAINS)


#get read lengths for all FASTQ files:
rule getReadLength:
	input: expand(config["FQPATH"] + "{{techname}}_{{capDesign}}_{{sizeFrac}}.fastq", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
#		config["FQPATH"] + "{techname}_{capDesign}_{sizeFrac}.fastq"
	output: "{techname}_{capDesign}_{sizeFrac}.readlength.tsv"
	shell: "fastq2fasta.pl {input} | FastaToTbl | awk -v s={wildcards.sizeFrac} '{{print s\"\\t\"length($2)}}' > {output}"

#aggregate read length data over all fractions of a given capDesign:
rule aggReadLength:
	input: expand("{{techname}}_{{capDesign}}_{sizeFrac}.readlength.tsv", filtered_product, techname=TECHNAMES, capDesign=CAPDESIGNS, sizeFrac=SIZEFRACS)
	#input: glob.glob(os.path.join("{capDesign}_*.readlength.tsv"))
	output: config["STATSDATADIR"] + "{techname}_{capDesign}_all.readlength.tsv"
	shell: "cat {input} > {output}"
