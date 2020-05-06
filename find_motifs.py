import argparse
parser = argparse.ArgumentParser()

from Bio import SeqIO
from Bio import Seq

parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
args = parser.parse_args()

def find_repeats(s, max_length, min_length):
   for i in range(len(s)):
      for j in range(min_length, max_length+1):
         count = 1
         while s[i:i+j] == s[i+j*count:i+j*count+j]: count += 1
         if count > 1:
               yield s[i:i+j], i, count

for record in SeqIO.parse(args.fasta, "fasta"):
      seqid = record.name 
      s = record.seq.upper()
      for pattern, position, count in find_repeats(s, 5, 3):
         print("%s \t %6s \t (%d, %d) \t %d" % (seqid, pattern, position, position + count*len(pattern), count))





