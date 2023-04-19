process EXT_TR{
    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-genomicalignments bioconda::r-dplyr bioconda::r-tibble conda-forge::r-tidyr bioconda::r-ggplot2 bioconda::r-viridis bioconda::grid"
    publishDir "$params.outdir/bamslam", mode: 'copy'

    input:
    path bam
    path gtf


    output:
    path "extended_transcriptome.fa", emit: extended_transcriptome

    script:
    """
    awk '\$3=="exon"' ${extended_annotation} > gtf_exon.gtf

    gffread -w extended_transcriptome.fa -g ${fasta} gtf_exon.gtf
    """
}