# nxf-minimapper
Nextflow pipeline to map long reads to a reference. Mainly for use in plasmid or amplicon verification.
The pipeline accepts a reference file (many file formats are supported) and a file or a directory of files
with sequencing reads. It performs alignment of the reads to the reference using `minimap2`, and 
generates an Integrated Genomic Viewer html report. The resulting alignments (`bam` files) can be also visualized in 
standalone IGV or other software.

## Running the pipeline
The pipeline can be run from the command line
```bash
nextflow run angelovangel/nxf-minimapper --ref --ref_format --fastq
```
or imported in EPI2ME and run there. For this, install EPI2ME 

## Output


