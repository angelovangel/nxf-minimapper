#!/usr/bin/env python3

# Convert various formats to fasta

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


in_file = sys.argv[1]
in_format = sys.argv[2]
out_file = sys.argv[3]

basename = os.path.basename(in_file).split('.', 1)[0]
n = 0

# do not rely on valid accession/header, take file name
# SeqIO.convert(in_file, in_format, out_file, "fasta")

for record in SeqIO.parse(in_file, in_format):
    n += 1
    id = record.id
    seq = record.seq
    record = SeqRecord(seq, id = basename + "_" + str(n), description = "")
    SeqIO.write(record, out_file, "fasta")
