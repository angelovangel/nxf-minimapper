#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO

##### parse arguments
parser = argparse.ArgumentParser(description="Converters a .fastq file to a .phd file.");
parser.add_argument('-i', '--input', required=True, help="The name of the input .fastq file.");
parser.add_argument('-o', '--output', required=True, help="The name of the file to output to (will overwrite file if it already exists)");

args = parser.parse_args();
#####

##### ___main___ #####

records = SeqIO.parse(args.input, "fastq")
count = SeqIO.write(records, args.output, "phd")
print("Converted %i records" % count)