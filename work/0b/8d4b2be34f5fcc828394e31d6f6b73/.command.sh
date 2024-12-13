#!/bin/bash -ue
# construct bed file
len=$(faster2 -l sample22.validated.fasta)
header=$(grep ">" sample22.validated.fasta | cut -c 2-)
echo -e "$header	0	$len	Primary alignments: X" > bedfile.bed


create_report         bedfile.bed         --fasta sample22.validated.fasta         --tracks sample22.bam         --output sample22.igvreport.html         --flanking 200
