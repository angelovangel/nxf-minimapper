#!/usr/bin/env python3

import sys
import argparse

##### parse arguments
parser = argparse.ArgumentParser(description="Converters a .vcf file to a .bed file.");
parser.add_argument('-i', '--input', required=True, help="The name of the input .vcf file.");
parser.add_argument('-o', '--output', required=False, help="The name of the file to output to (will overwrite file if it already exists). If this argument is not specified, output will be to the input file name with a .bed extension");

_args = parser.parse_args();
#####

##### ___main___ #####

if _args.input[-4:] != ".vcf":
    print("Input file is not a .vcf file");
    sys.exit();

if _args.input == _args.output:
    print("Input and output cannot be the same file")
    sys.exit()

vcf_file = open(_args.input)
if _args.output != None:
    bed_file_name = _args.output
else:
    bed_file_name = _args.input[:-4] + ".bed"
bed_file = open(bed_file_name, "w");

for line in vcf_file:
    # Skip header lines
    if line.startswith('#'):
        continue
        
    try:
        firstTabIndex = line.index('\t')
        secondTabIndex = line.index('\t', firstTabIndex + 1)
        pos = int(line[firstTabIndex+1:secondTabIndex])
        rest = line[secondTabIndex + 1:].replace('\t', '/')
        print(rest)
        bed_file.write('\t'.join([line[0:firstTabIndex], str(pos - 1), str(pos), rest]))
        #bed_file.write('\t'.join([line[0:firstTabIndex], str(pos - 1), str(pos), line[secondTabIndex + 1:]]))
    except ValueError as e:
        print(f"Error processing line: {line.strip()}")
        print("Line appears to be malformed. Please ensure it follows VCF format.")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error processing line: {line.strip()}")
        print(f"Error: {str(e)}")
        sys.exit(1)

print("Outputed to " + bed_file_name);

#########################