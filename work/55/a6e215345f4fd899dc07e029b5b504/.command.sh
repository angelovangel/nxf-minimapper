#!/bin/bash -ue
minimap2 -a -x lr:hq         sample22.validated.fasta         sample22.fastq.gz         | samtools view -S -b -         | samtools sort -o sample22.bam -
samtools index sample22.bam
