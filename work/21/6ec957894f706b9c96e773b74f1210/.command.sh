#!/bin/bash -ue
create_report         bedfile.bed         --fasta sample22.validated.fasta         --tracks sample22.bam         --output sample22.igvreport.html         --flanking 200
