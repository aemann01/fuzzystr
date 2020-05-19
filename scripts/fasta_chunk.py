#!/usr/bin/env python3

'''
Use: python3 fasta_chunk.py input.fa 
'''

import sys
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = False
            if not entry:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

record_iter = SeqIO.parse(sys.argv[1], 'fasta')

for i, batch in enumerate(batch_iterator(record_iter, int(sys.argv[2])), start=1):
    filename = 'group_{}.fasta'.format(i)
    count = SeqIO.write(batch, filename, 'fasta')
    print('Wrote {} records to {}'.format(count, filename))