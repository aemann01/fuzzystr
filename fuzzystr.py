#!/usr/bin/env python3

'''
Fuzzy match primer protocol for STR project. Use: python3 fuzzystr.py -f input.fasta -p primerFile.txt -o output.txt
'''

import pandas as pd
import numpy as np
import os
import argparse
parser = argparse.ArgumentParser()

from fuzzysearch import find_near_matches
from Bio import SeqIO

parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-p', '--primers', help='list of primers', required=True)
parser.add_argument('-o', '--output', help='output file', default='output.txt')
parser.add_argument('-e', '--edit', help='desired edit distance cutoff (float)', default='0.25')
args = parser.parse_args()

def fuzzystr():
	# initialize empty dataframe
	df = pd.DataFrame(columns=["seqid", "seq", "trim", "motif", "length", "trimlen", "locus", "fwstart", "fwend", "rvstart", "rvend", "fwdist", "rvdist", "fwmatch_seq", "rvmatch_seq"])
	record_count = 0
	for record in SeqIO.parse(args.fasta, "fasta"):
		seqid = record.name 
		seq = record.seq
		length = len(record.seq)
		record_count += 1

		# print record count to screen every 100 loops
		if record_count % 100 == 0:
			print(record_count)

		# iterate over primer file for each record
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

				# break if either forward or reverse match are not found
				if len(fwmatch) == 0 or len(rvmatch) == 0:
					break
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
					break
				# trim region
				else:
					rvend += 1
					trim = seq[fwstart:rvend]
				trimlen = len(trim)
					
				#check if motif is found at least once in the sequence
				ed3 = round(len(motif)*float(args.edit))
				motif_check = find_near_matches(motif, trim, max_l_dist=ed3)
				if motif_check == 0:
					break
			# populate dataframe
			df = df.append({"seqid": seqid, "seq": seq, "trim": trim, "trimlen": trimlen, "motif": motif, "length": length, "locus": locus, "fwstart": fwstart, "fwend": fwend, "rvstart": rvstart, "rvend": rvend, "fwdist": fwdist, "rvdist": rvdist, "fwmatch_seq": fwseq, "rvmatch_seq": rvseq}, ignore_index=True)
	# write to file
	with open(args.output, "w") as output:
		df.to_csv(output, index=False)

fuzzystr()



























