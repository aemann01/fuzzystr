#!/usr/bin/env python3

'''
Count observed motifs after trimming using fuzzystr
'''




def motifcount():
	res = [seq[i:j] for i in range(len(seq)) for j in range(i+1, len(seq) + 1)]
	countdict = {i:res.count(i) for i in res