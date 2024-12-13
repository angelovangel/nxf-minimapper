if (params.help) {
    helpMessage()
    exit(0)
}



ref_ch = Channel.fromPath(params.ref, type: 'file', checkIfExists: true)
reads_ch = Channel.fromPath(params.fastq, type: 'file', checkIfExists: true)

// evtl convert .dna to fasta, check for valid ref
process VALIDATE_REF {
    container 'aangeloo/nxf-tgs:latest'

    input: path(ref)
    output: path("${ref.simpleName}.validated.fasta"), emit: validated_ref_ch

    script:
    """
    seqkit seq -v $ref > ${ref.simpleName}.validated.fasta
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
    # construct bed file
    len=\$(faster2 -l ${ref})
    header=\$(grep ">" ${ref} | cut -c 2-)
    echo -e "\$header\t0\t\$len\tPrimary alignments: X" > bedfile.bed


    create_report \
        bedfile.bed \
        --fasta ${ref} \
        --tracks ${bam} \
        --output ${bam.simpleName}.igvreport.html \
        --flanking 200 \
    """
}

workflow {
    VALIDATE_REF(ref_ch)
    MINIMAP(VALIDATE_REF.out.validated_ref_ch, reads_ch)
    
    VALIDATE_REF.out.validated_ref_ch \
    .combine(MINIMAP.out.bam_ch) \
    | IGV
}