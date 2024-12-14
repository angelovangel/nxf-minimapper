

ref_ch = Channel.fromPath(params.ref, type: 'file', checkIfExists: true)
reads_ch = Channel.fromPath(params.fastq, type: 'any', checkIfExists: true)

// evtl convert .dna to fasta, check for valid ref
process VALIDATE_REF {
    //container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input: path(ref)
    output: path("${ref.simpleName}.validated.fasta"), emit: validated_ref_ch

    script:
    """
        convert2fasta.py $ref "${params.format}" ${ref.simpleName}.validated.fasta
    """
}

// // do actual mapping
process MINIMAP {
    container 'aangeloo/nxf-tgs:latest'
    publishDir "$params.outdir", mode: 'copy'

    input:
        path(ref)
        path(fastq)
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

// generate IGV report
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
    echo -e "\$header\t0\t\$len\tPrimary alignments: \$percent_aln of \$all_reads reads" > bedfile.bed


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
    MINIMAP(VALIDATE_REF.out.validated_ref_ch, reads_ch)
    
    VALIDATE_REF.out.validated_ref_ch \
    .combine(MINIMAP.out.bam_ch) \
    | IGV
}