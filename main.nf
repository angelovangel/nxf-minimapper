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
    ref           : path to a reference file, containing one sequence record
    format        : file format of reference (one of 'fasta', 'genbank', 'embl', 'snapgene')
    minimap_args  : additional command-line arguments passed to minimap2, pass these as strings '--param1 value --param2 value'...
    outdir        : where to save results, default is 'output'
    variants      : whether to include variants in the IGV reports, default is false
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
    variants        : ${params.variants}
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
    publishDir "$params.outdir/00-alignments/", mode: 'copy'

    input: path(ref)
    output: path("${ref.simpleName}.reference.fasta"), emit: validated_ref_ch

    script:
    """
    convert2fasta.py $ref "${params.format}" ${ref.simpleName}.reference.fasta
    """
}


process MINIMAP {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir/00-alignments/", mode: 'copy'

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

process MEDAKA_VARIANT {
    container 'ontresearch/medaka:latest'
    publishDir "$params.outdir/02-variants/", mode: 'copy'

    input: tuple path(ref), path(fastq)
    output: path("*variants.{vcf,vcf.gz}"), emit: vcf_ch

    script:
    """
    medaka_variant -i $fastq -r $ref 
    mv medaka/medaka.annotated.vcf ${fastq.simpleName}.variants.vcf
    """
}

process MEDAKA_CONSENSUS {
    container 'ontresearch/medaka:latest'
    publishDir "$params.outdir/00-alignments/", mode: 'copy'

    input: tuple path(ref), path(fastq)
    output: tuple val(sample), path("*consensus.fastq"), emit: consensus_ch

    script:
    sample = fastq.simpleName
    """
    # fill gaps with N in consensus
    medaka_consensus -r N -q -i $fastq -d $ref 
    mv medaka/consensus.fastq ${fastq.simpleName}.consensus.fastq
    """
}

// get coverage statistics for all samples that were mapped to ref
process COVERAGE_STATS {
    container 'aangeloo/nxf-tgs:latest'
    //publishDir "$params.outdir", mode: 'copy'

    input: path(bam)
    output: tuple val(sample), path("coverage.tsv"), emit: coverage_ch

    script:
    sample = bam.simpleName
    
    // record samtools depth as list col with "|" sep2 in coverage.tsv
    """
    # determine sampling nth so that approx 300 points are collected for sparkline depth of coverage
    endpos=\$(samtools coverage -H ${bam} | cut -f3)
    nth=\$(awk -v var="\$endpos" 'BEGIN {print int(var / 300) }')
    echo "nth: \$nth"
    
    samtools coverage -H ${bam} | awk '{print "$sample\t" \$0}' - > coverage
    samtools depth -aa ${bam} | awk -v var="\$nth" 'NR % var == 0' - | cut -f2,3 | tr '\n' '|' | tr '\t' ':' > depth
    paste -d "\t" coverage depth > coverage.tsv
    """
}

process CONSENSUS_STATS {
    container 'pegi3s/biopython:latest'

    input: tuple val(sample), path(consensus), path(coverage)
    output: path("coverage2.tsv"), emit: spark_ch

    script:
    """
    # convert consensus.fastq to phd and add to coverage.tsv
    fastq2phd.py -i $consensus -o ${consensus.simpleName}.consensus.phd
    cat ${consensus.simpleName}.consensus.phd | grep -v -e "BEGIN_*" -e "END*" > phd.file
    
    # determine sampling nth so that approx 300 points are collected for sparkline depth of coverage
    endpos=\$(cat phd.file | wc -l | awk '\$1=\$1' -)
    nth=\$(awk -v var="\$endpos" 'BEGIN {print int(var / 300) }')
    # average cons q-score over 5 bases
    run-average.awk col=2 size=10 phd.file | awk -v var="\$nth" 'NR % var == 0' - | tr '\n' '|' > cons
    paste -d "\t" $coverage cons > coverage2.tsv
    """
}

process COVERAGE_SUMMARY {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input:  path("coverage2.tsv")
    output: path("00-alignment-summary*")

    script:
    """
    coverage_summary.R "*.tsv" $workflow.runName
    """

}


process IGV {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir/01-igv-reports/", mode: 'copy'
    
    input: 
    tuple val(sample), path(bam), path(vcf), path(ref)
    
    output:
    path "*.igvreport.html"

    script:
    def vcf_track = params.variants ? "$vcf" : ""
    """
    # dynamic calculation for subsampling, subsample for > 500 alignments
    count=\$(samtools view -c ${bam[0]})
    percent_aln=\$(samtools flagstat ${bam[0]} | grep 'primary mapped' | cut -d"(" -f 2 | cut -d" " -f1)
    all_reads=\$(samtools flagstat ${bam[0]} | grep 'primary\$' | cut -d" " -f1)
    subsample=\$(echo \$count | awk '{if (\$1 <500) {print 1} else {print 500/\$1}}')

    # construct bed file
    len=\$(faster2 -l ${ref})
    header=\$(grep ">" ${ref} | cut -c 2-)
    echo -e "\$header\t0\t\$len\tShowing \$subsample fraction of \$all_reads alignments" > bedfile.bed

    # add vcf to bed file
    if [ -n "$vcf_track" ]; then
        vcf2bed.py -i $vcf -o vcf.bed
        head -n 10 vcf.bed >> bedfile.bed
    fi

    create_report \
        bedfile.bed \
        --fasta ${ref} \
        --tracks ${vcf_track} ${bam[0]} \
        --output ${sample}.igvreport.html \
        --flanking 200 \
        --subsample \$subsample
    """
}

workflow {
    VALIDATE_REF(ref_ch)

    VALIDATE_REF.out.validated_ref_ch
    .combine(reads_ch) \
    | (MINIMAP & MEDAKA_VARIANT & MEDAKA_CONSENSUS) 

    MINIMAP.out.bam_ch
    .flatten()
    .filter( ~/.*bam$/ )
    //.view()
    | COVERAGE_STATS

    MEDAKA_CONSENSUS.out.consensus_ch
    .join(COVERAGE_STATS.out.coverage_ch, by:0)
    .set { stats_ch }

    stats_ch
    | CONSENSUS_STATS
    
    CONSENSUS_STATS.out.spark_ch
    .collect()
    //.view()
    | COVERAGE_SUMMARY

    MINIMAP.out.bam_ch \
    .map{ [ it.toString().split("/").last().split("\\.")[0], it ] } //make keys for join
    .set { minimap_ch }
    
    MEDAKA_VARIANT.out.vcf_ch \
    .map{ [ it.toString().split("/").last().split("\\.")[0], it ] } //make keys for join
    .set { medaka_ch }

    minimap_ch \
    .join(medaka_ch, by: 0) \
    .combine(VALIDATE_REF.out.validated_ref_ch) // validated_ref is always 1 file
    //.view()
    | IGV

}