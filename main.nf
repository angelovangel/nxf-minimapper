#!/usr/bin/env nextflow

if (params.help) {
    helpMessage()
    exit(0)
}

def helpMessage() {
    log.info """\
    ========================================================================================================================
    NXF - MINIMAPPER
    Map long reads to reference (for plasmid/amplicon verification)
    ========================================================================================================================
    Parameters:
    -----------------------------------
    fastq         : path to a fastq file or a directory of fastq files to be mapped to the reference
    ref           : path to a reference file
    format        : file format of reference (one of 'fasta', 'genbank', 'embl', 'snapgene')
    minimap_args  : additional command-line arguments passed to minimap2, pass these as strings '--param1 value --param2 value'...
    outdir        : where to save results, default is 'output'
    """
    .stripIndent(true)
}

log.info """\
    ========================================================================================================================
    NXF - MINIMAPPER
    Map long reads to reference (for plasmid/amplicon verification)
    ========================================================================================================================
    fastq           : ${params.fastq}
    ref             : ${params.ref}
    format          : ${params.format}
    outdir          : ${params.outdir}
    minimap_args    : ${params.minimap_args}
    """
    .stripIndent(true)

ref_ch = Channel.fromPath(params.ref, type: 'file', checkIfExists: true)

// want to have either a single fastq file or a directory with fastq files as a query
// this is "path" in the schema to allow selecting both file and dir in EPI2ME
// the channel emits files only
myfile = file(params.fastq)
mypattern = "*.{fastq,fastq.gz,fq,fq.gz}"

if ( myfile.isDirectory() ) {
    reads_ch = Channel.fromPath(params.fastq + "/" + mypattern, type: 'file', checkIfExists: true)
} else {
    reads_ch = Channel.fromPath(params.fastq, type: 'file', checkIfExists: true)
}

process VALIDATE_REF {
    container 'pegi3s/biopython:latest'
    publishDir "$params.outdir", mode: 'copy'

    input: path(ref)
    output: path("${ref.simpleName}.reference.fasta"), emit: validated_ref_ch

    script:
    """
    convert2fasta.py $ref "${params.format}" ${ref.simpleName}.reference.fasta
    """
}


process MINIMAP {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input: tuple path(ref), path(fastq)
    output: path("*.{bam,bai}"), emit: bam_ch

    script:
    """
    minimap2 -a -x lr:hq \
        $ref \
        $fastq \
        | samtools view -S -b - \
        | samtools sort -o ${fastq.simpleName}.bam -
    samtools index ${fastq.simpleName}.bam
    """
}

// get coverage statistics for all samples that were mapped to ref
process COVERAGE_STATS {
    container 'aangeloo/nxf-tgs:latest'
    //publishDir "$params.outdir", mode: 'copy'

    input: path(bam)
    output: path("coverage.tsv"), emit: coverage_ch

    script:
    sample = bam.simpleName
    """
    samtools coverage -H ${bam} | awk '{print "$sample\t" \$0}' - > coverage.tsv
    """
}

process COVERAGE_SUMMARY {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input:  path("coverage*.tsv")
    output: path("00-alignment-summary*")

    script:
    """
    coverage_summary.R "*.tsv" $workflow.runName
    """

}


process IGV {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'
    
    input: 
    tuple path(ref), path(bam), path(bai)
    
    output:
    path "*.igvreport.html"

    script:
    """
    # dynamic calculation for subsampling, subsample for > 500 alignments
    count=\$(samtools view -c ${bam})
    percent_aln=\$(samtools flagstat ${bam} | grep 'primary mapped' | cut -d"(" -f 2 | cut -d" " -f1)
    all_reads=\$(samtools flagstat ${bam} | grep 'primary\$' | cut -d" " -f1)
    subsample=\$(echo \$count | awk '{if (\$1 <500) {print 1} else {print 500/\$1}}')

    # construct bed file
    len=\$(faster2 -l ${ref})
    header=\$(grep ">" ${ref} | cut -c 2-)
    echo -e "\$header\t0\t\$len\tShowing \$subsample fraction of \$all_reads alignments" > bedfile.bed


    create_report \
        bedfile.bed \
        --fasta ${ref} \
        --tracks ${bam} \
        --output ${bam.simpleName}.igvreport.html \
        --flanking 200 \
        --subsample \$subsample
    """
}

workflow {
    VALIDATE_REF(ref_ch)
    
    VALIDATE_REF.out.validated_ref_ch \
    .combine(reads_ch)
    | MINIMAP

    VALIDATE_REF.out.validated_ref_ch \
    .combine(MINIMAP.out.bam_ch) \
    | IGV

    MINIMAP.out.bam_ch
    .flatten()
    .filter( ~/.*bam$/ )
    //.view()
    | COVERAGE_STATS
    
    COVERAGE_STATS.out.coverage_ch
    .collect()
    //.view()
    | COVERAGE_SUMMARY
    
}