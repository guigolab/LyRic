#!/usr/bin/env python


# Modified from: <<<
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# Modified from >>>

# Uses pyfaidx to extract FASTA subsequences (Shirley MD, Ma Z, Pedersen B, Wheelan S. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196. 2015.)

"""
ScoreSpliceCandidate.py

Given the context of a GT or AG in a string of DNA bases and precomputed coefficients
files, return a score that indicates how likely it is to be a splice donor or acceptor,
respectively, using the maximum entropy method [1].

Usage: create a DonorPredictor or AcceptorPredictor class instance using a coefficients
    file, and then invoke one or more times supplying 3 bases before the GT and 4 bases
    after for donor predictions or 18 bases before the AG and 3 bases after for acceptor
    predictions.

Example:
    from ScoreSpliceCandidate import DonorPredictor, AcceptorPredictor
    donorPredictor = DonorPredictor('Hsap.donor.mecoef')
    print(donorPredictor('CAG', 'GAGC')) # Score for CAG-GT-GAGC
    # Prints 10.173901761046155
    acceptorPredictor = AcceptorPredictor('Hsap.acceptor.mecoef')
    print(acceptorPredictor('TGTGTGCCTTTCACTTTC', 'GCT')) # Score for TGTGTGCCTTTCACTTTC-AG-GCT
    # Prints 10.791433170214455

See comment at the bottom of this file for format of coefficient files.

[1] Yeo, G., & Burge, C. B. (2004). Maximum entropy modeling of short sequence motifs with
applications to RNA splicing signals. Journal of Computational Biology : a Journal of
Computational Molecular Cell Biology, 11(2-3), 377-394. doi:10.1089/1066527041410418
"""
from __future__ import division, print_function
import os, struct, math, sys
from Bio import SeqIO
from pyfaidx import Fasta

assert sys.version_info[0] == (2)


class DonorPredictor(object) :
    def __init__(self, donorCoefFileName) :
        self.file = open(os.path.abspath(os.path.expanduser(donorCoefFileName)))
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 16384, \
            '%s is not a valid donor coefficients file.' % donorCoefFileName

    def __call__(self, prev3bases, next4bases) :
        # Score potential GT splice donor site given 3 prev bases and 4 next ones.
        assert len(prev3bases) == 3, len(prev3bases)
        assert len(next4bases) == 4, len(next4bases)
        index = _bases_to_number(prev3bases + next4bases)
        self.file.seek(index * RecordSize)
        coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
        return math.log(16.302010666666664 * coeff, 2)

class AcceptorPredictor(object) :
    def __init__(self, acceptorCoefFileName) :
        self.file = open(os.path.abspath(os.path.expanduser(acceptorCoefFileName)))
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 82560, \
            '%s is not a valid acceptor coefficients file.' % acceptorCoefFileName

    def __call__(self, prev18bases, next3bases) :
        # Score potential AG splice acceptor site given 18 prev bases and 3 next ones.
        assert len(prev18bases) == 18, len(prev18bases)
        assert len(next3bases) == 3, len(next3bases)
        bases = prev18bases + next3bases
        coeffsCombination = 1
        for ii, (start, end) in enumerate(AcceptorStartEnds) :
            index = _bases_to_number(bases[start : end + 1])
            self.file.seek((AcceptorArrayLengthSums[ii] + index) * RecordSize)
            coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
            if ii < 5 :
                coeffsCombination *= coeff
            else :
                coeffsCombination /= coeff
        return math.log(16.3482025 * coeffsCombination, 2)

def _bases_to_number(bases) :
    "Convert a string of DNA bases to a base-4 number."
    BaseMap = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    result = 0
    for b in bases :
        result = 4 * result + BaseMap[b]
    return result

RecordFormat = '<d' # little-endian double
RecordSize = struct.calcsize(RecordFormat) # Typically 8

# Relevant intervals for acceptor site prediction
AcceptorStartEnds = [(0, 6), (7, 13), (14, 20), (4, 10), (11, 17),
                     (4, 6), (7, 10), (11, 13), (14, 17)]
# Lengths of the coefficient arrays for acceptor sites.
AcceptorArrayLengths = [4 ** (end - start + 1) for start, end in AcceptorStartEnds]
AcceptorArrayLengthSums = [sum(AcceptorArrayLengths[:ii])
                           for ii in range(len(AcceptorArrayLengths))]

gffFile = sys.argv[1]
genomeFa= sys.argv[2]
donorCoefFileName = sys.argv[3]
acceptorCoefFileName = sys.argv[4]
donorPredictor = DonorPredictor(donorCoefFileName)
acceptorPredictor = AcceptorPredictor(acceptorCoefFileName)
chromosomes = Fasta(genomeFa)#, read_ahead=10000)

minDonorScore=-19
minAcceptorScore=-24
if gffFile == "-":
    gff = sys.stdin
else:
    gff = open(gffFile)

gffdata=gff.readlines()
line_number=0
for line in gffdata:
    line_number+=1
    line=line.rstrip()
    fields= line.split("\t")
    if len(fields) != 9:
        print("Invalid number of fields in " + gffFile + " line number " + str(line_number) + ", must be 9. Is this GFF/GTF?\n", file=sys.stderr)
        sys.exit(1)
    if fields[2] == 'intron':
        chrom=fields[0]
        leftStart=int(fields[3])-1
        rightStart=int(fields[4])-2
        strand=fields[6]
        leftEnd=leftStart+2
        rightEnd=rightStart+2
        if strand == '+' or strand == '-':
            if strand == '+':
                donor = chromosomes[chrom][leftStart-3:leftEnd+4].seq.upper()
                acceptor = chromosomes[chrom][rightStart-18:rightEnd+3].seq.upper()
            elif strand == '-':
                donor = chromosomes[chrom][rightStart-4:rightEnd+3].reverse.complement.seq.upper()
                acceptor = chromosomes[chrom][leftStart-3:leftEnd+18].reverse.complement.seq.upper()
            donorDiNt= donor[3:5]
            acceptorDiNt=acceptor[18:20]
            if (donorDiNt == 'GT' and acceptorDiNt == 'AG' and len(donor) == 9 and len(acceptor) == 23 and 'N' not in donor and 'N' not in acceptor):
                donorFlank1=donor[:3]
                donorFlank2=donor[5:]
                acceptorFlank1=acceptor[:18]
                acceptorFlank2=acceptor[20:]
                #print(donorFlank1 + ' ' + donorFlank2)
                #print(acceptorFlank1 + ' ' + acceptorFlank2)
                scoreDonor=donorPredictor(donorFlank1, donorFlank2)
                scoreAcceptor=acceptorPredictor(acceptorFlank1, acceptorFlank2)
            else:
                scoreAcceptor = minAcceptorScore
                scoreDonor = minDonorScore
        else:
            scoreAcceptor = minAcceptorScore
            scoreDonor = minDonorScore
        fields[8] = fields[8] + " accMeScore \"" + str(scoreAcceptor) + "\"; donMeScore \"" + str(scoreDonor) + "\";"
        print ("\t".join(fields))
        #print (str(scoreAcceptor) + " " + str(scoreDonor))

#print(chromosomes['chr13'][57632919:57633018].seq)


"""
Format of coefficient files:

The donor coefficient file consists of 4^7 numbers stored in little-endian double format.
The index into this file for a particular splice candidate is computed from the string
of 3 bases 5' of the GT concatenated to the 4 bases 3' of the GT. The index is a base-4
number created from the indices of the bases of this base string in ACGT with the first
being most significant. For example, the index for a potential splice site context of:
ACG-GT-TACG is 0*4^6 + 1*4^5 + 2*4^4 + 3*4^3 + 0*4^2 + 1*4^1 + 2*4^0
Once the coefficient is extracted, the score is
	log(16.302010666666664 * coeff, 2)

The acceptor coefficient file consists of 9 concatenated sequences containing
4^7, 4^7, 4^7, 4^7, 4^7, 4^3, 4^4, 4^3, and 4^4 numbers. The relevant context
is 18 bases 5' of the AG concatenated to 3 bases 3'. The 9 coefficients are extracted
from these arrays using indices computed from the substrings starting at the following
offsets with the following lengths:
0 7, 7 7, 14 7 (including the 3 bases 3' of the AG), 4 7, 11 7, 4 3, 7 4, 11 3, 14 4
Once these 9 coefficients are found, the score is:
    log(16.3482025 * (prod 1st 5 coeffs) / (prod other 4 coeffs), 2)
"""
