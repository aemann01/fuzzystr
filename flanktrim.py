 #!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()

from fuzzysearch import find_near_matches
from Bio import SeqIO

parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-fw', '--fwflank', help='forward flanking region to trim', required=True)
parser.add_argument('-rv', '--rvflank', help='reverse flanking region to trim', required=True)
parser.add_argument('-e', '--edit', help='desired edit distance cutoff (float)', default='0.25')
args = parser.parse_args()

def flanktrim():
	for record in SeqIO.parse(args.fasta, "fasta"):
		seqid = record.name 
		seq = record.seq
		fwflank = args.fwflank
		rvflank = args.rvflank
		ed1 = round(len(fwflank)*float(args.edit))
		ed2 = round(len(rvflank)*float(args.edit))
		fwmatch = find_near_matches(fwflank, seq, max_l_dist=ed1)
		rvmatch = find_near_matches(rvflank, seq, max_l_dist=ed2)
		
		# continue with next locus if either forward or reverse match are not found
		if len(fwmatch) == 0 or len(rvmatch) == 0:
			continue
		# for both forward and reverse, pick index of lowest edit distance match
		else:
			fwindex = []
			for i in range(0, len(fwmatch)):
				fwindex.append(fwmatch[i].dist)
			fwstart = fwmatch[fwindex.index(min(fwindex))].start
			fwend = fwmatch[fwindex.index(min(fwindex))].end
			fwdist = fwmatch[fwindex.index(min(fwindex))].dist
			fwseq = fwmatch[fwindex.index(min(fwindex))].matched
	
			rvindex = []
			for i in range(0, len(rvmatch)):
				rvindex.append(rvmatch[i].dist)
			rvstart = rvmatch[rvindex.index(min(rvindex))].start
			rvend = rvmatch[rvindex.index(min(rvindex))].end
			rvdist = rvmatch[rvindex.index(min(rvindex))].dist
			rvseq = rvmatch[rvindex.index(min(rvindex))].matched

		# break if reverse is downstream to forward primer
		if rvend < fwstart:
			continue
		# trim region
		else:
			trim = seq[fwstart:rvend]
			trimlen = len(trim)

		# check if trimmed sequence is within the fragment length range
		if trimlen <= 20 or trimlen >= 300:
			continue
		print(seqid, "\t", trim, "\t", trimlen, "\t", fwdist, "\t", rvdist, "\t")

def main():
	assert os.path.exists(args.fasta), 'Error! File does not exist: %s. Is the path correct?' % args.fasta
	flanktrim()

main()

__author__ = "Allison E. Mann"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "allison.e.mann@gmail.com"