process EXT_TR{
    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hd03093a_2' :
        'biocontainers/gffread:0.12.7--hd03093a_2' }"
    publishDir "$params.outdir/bamslam", mode: 'copy'

    input:
    path extended_annotation
    path fasta

    output:
    path "extended_transcriptome.fa", emit: extended_transcriptome
    path "versions.yml" , emit: versions
    
    script:
    """
    awk '\$3=="exon"' ${extended_annotation} > gtf_exon.gtf
    gffread -w extended_transcriptome.fa -g ${fasta} gtf_exon.gtf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}