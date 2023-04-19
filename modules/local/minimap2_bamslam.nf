process EXT_TR{
    tag "$meta.id"
    conda "bioconda::minimap2 bioconda::samtools"
    publishDir "$params.outdir/bamslam", mode: 'copy'

    input:
    path extended_transcriptome
    path fastq

    output:
    path "*.bam", emit: bam

    script:
    """
    minimap2 -t 16 -ax map-ont ${extended_transcriptome} ${fastq} > ${meta.id}.sam

    samtools view -bh -o ${meta.id}.bam ${meta.id}.sam
    """
}