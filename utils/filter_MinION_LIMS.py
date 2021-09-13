#!/usr/bin/env python

import sys
import pandas as pd

inputF = sys.argv[1]

if inputF == '-':
    sampleAnnot = pd.read_table(sys.stdin)
else:
    sampleAnnot = pd.read_table(inputF)


subsetSampleAnnot = sampleAnnot[sampleAnnot['filesystem_path'].str.match(
    '.*/.*') == True]  # select rows where column 'filesystem_path' looks like a Unix path (https://kanoki.org/2019/11/12/how-to-use-regex-in-pandas/)
subsetSampleAnnot = subsetSampleAnnot[[
    'Experiment_group', 'Flowcell_product_code', 'Kit_ID', 'filesystem_path', 'fast5_subdir']]

subsetSampleAnnot.to_csv(sys.stdout, sep='\t', header=False, index=False)
