# nxf-minimapper
Nextflow pipeline to map long reads to a reference. Mainly for use in plasmid or amplicon verification.
The pipeline accepts a reference file containing one record and either:
- a fastq file or 
- a directory of fastq files
  
It is assumed that each sequence file cotains the reads from one sample that are to be mapped to the reference. It performs alignment of the reads to the reference using `minimap2`, and 
generates an Integrated Genomic Viewer html report. The resulting alignments (`bam` files) can be also visualized in 
standalone IGV or other software.  
Supported reference file formats:
 - fasta
 - genbank
 - EMBL
 - Snapgene

## Running the pipeline
Only `nextflow ` and `docker` are needed. The pipeline can be run from the command line
```bash
nextflow run angelovangel/nxf-minimapper --ref reference.dna --format snapgene --fastq sample1.fastq
```

The pipeline can also be imported in EPI2ME and run there. For this, [install EPI2ME](https://labs.epi2me.io/epi2me-docs/quickstart/) and import the pipeline using this link:  
`https://github.com/angelovangel/nxf-minimapper` 

## Output
Results are saved in a folder named `output`. A `bam` file is generated for every sample (in `output/00-alignments`) that can be opened in IGV. Also, an alignment summary and IGV report html files are generated. To open the alignments in IGV:
 - install IGV for your platform (https://igv.org/doc/desktop/#DownloadPage/) 
 - open reference file - `Genomes` - `Load genome from File` - select the reference.fasta file from `output/00-alignments`
 - open alignment - `File` - `Load from File` - select one or more `bam` files
  
If you use EPI2ME, all the reports are visible directly there under `Reports`.

