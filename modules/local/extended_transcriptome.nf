process EXT_TR{
    conda "bioconda::gffread=0.12.7"
    publishDir "$params.outdir/bamslam", mode: 'copy'

    input:
    path extended_annotation
    path fasta

    output:
    path "extended_transcriptome.fa", emit: extended_transcriptome

    script:
    """
    awk '\$3=="exon"' ${extended_annotation} > gtf_exon.gtf

    gffread -w extended_transcriptome.fa -g ${fasta} gtf_exon.gtf
    """
}