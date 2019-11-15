import glob
from collections import defaultdict
import os
import itertools
import sys

# do not rely on FASTQ filenames, otherwise we need to rename them each time we add some metadata

rule all:
    input:
        expand({library}.bam, )




