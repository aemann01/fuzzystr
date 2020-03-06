#!/usr/bin/env python3

'''
Fuzzy match primer protocol for STR project. Use: python3 test_script.py -f input.fastq -p primerFile.txt -o output.txt
'''

import pandas as pd
import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()

from fuzzysearch import find_near_matches
from Bio import SeqIO

parser.add_argument('-f', '--fastq', help='input fastq file', required=True)
parser.add_argument('-p', '--primers', help='list of primers', required=True)
parser.add_argument('-o', '--output', help='output file', default='output.txt')
parser.add_argument('-e', '--edit', help='desired edit distance cutoff (float)', default='0.25')
args = parser.parse_args()

def fuzzystr():
	df = pd.DataFrame(columns=["seqid", "seq", "trim", "motif", "length", "trimlen", "locus", "fwstart", "fwend", "rvstart", "rvend", "fwdist", "rvdist", "fwmatch_seq", "rvmatch_seq"])
	record_count = 0
	for record in SeqIO.parse(args.fastq, "fastq"):
		seqid = record.name 
		seq = record.seq
		length = len(record.seq)
		record_count += 1

		if record_count % 100 == 0:
			print(record_count)

		with open(args.primers, "r") as f:
			for line in f:
				fields = line.strip().split()
				locus = fields[0]
				fwprimer = fields[2]
				rvprimer = fields[3]
				motif = fields[4]
				ed1 = round(len(fwprimer)*float(args.edit))
				ed2 = round(len(rvprimer)*float(args.edit))
				fwmatch = find_near_matches(fwprimer, seq, max_l_dist=ed1)
				rvmatch = find_near_matches(rvprimer, seq, max_l_dist=ed2)

				if fwmatch == []:
					break
				else:
					distances = []
					for i in range(0, len(fwmatch)):
						distances.append(fwmatch[i].dist)
					fwstart = fwmatch[distances.index(min(distances))].start
					fwend = fwmatch[distances.index(min(distances))].end
					fwdist = fwmatch[distances.index(min(distances))].dist
					fwseq = fwmatch[distances.index(min(distances))].matched

				if rvmatch == []:
					rvhit = "NA"
					rvstart = "NA"
					rvend = "NA"
					rvdist = "NA"
					rvseq = "NA"
				else:
					distances = []
					for i in range(0, len(rvmatch)):	
						distances.append(rvmatch[i].dist)
					rvstart = rvmatch[distances.index(min(distances))].start
					rvend = rvmatch[distances.index(min(distances))].end
					rvdist = rvmatch[distances.index(min(distances))].dist
					rvseq = rvmatch[distances.index(min(distances))].matched

				if rvstart == "NA":
					trim = seq[fwstart:]
				elif fwstart < rvend:
					rvend += 1
					trim = seq[fwstart:rvend]
				else:
					fwend += 1
					trim = seq[rvstart:fwend]
				trimlen = len(trim)

			df = df.append({"seqid": seqid, "seq": seq, "trim": trim, "motif": motif, "length": length, "trimlen": trimlen, "locus": locus, "fwstart": fwstart, "fwend": fwend, "rvstart": rvstart, "rvend": rvend, "fwdist": fwdist, "rvdist": rvdist, "fwmatch_seq": fwseq, "rvmatch_seq": rvseq}, ignore_index=True)

	with open(args.output, "w") as output:
		df.to_csv(output, index=False)

fuzzystr()

#def motifcount():
	#res = [seq[i:j] for i in range(len(seq)) for j in range(i+1, len(seq) + 1)]
	#countdict = {i:res.count(i) for i in res




#def alignref():

#def filter():
## TO DO: length and motif filtering function

#def strstats():
## TO DO: get the rep # and motif information, generate statistics


















