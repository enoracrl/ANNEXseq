process ISOFORMSWITCHANALYZER {

    conda "bioconda::bioconductor-isoformswitchanalyzer bioconda::bioconductor-data.table bioconda::r-genomicfeatures conda-forge::r-dplyr=1.0.10 conda-forge::r-tidyr conda-forge::r-ggvenn conda-forge::r-cowplot conda-forge::r-upsetr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.0.8--r42hc247a5b_0' :
        'quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
    publishDir "$params.outdir/IsoformSwitchAnalyzeR", mode: 'copy'

    input:
        file matrix
        file design
        file comparisons
        path transcriptome
        path gtf

    script:
    """
    IsoformSwitchAnalyzeR.r \\
        --annotation=$gtf \\
        --transcriptome=$transcriptome \\
        --design=$design \\
        --count_matrix=$matrix \\
        --comparisons=$comparisons

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-isoformswitchanalyzer: \$(Rscript -e "library(isoformswitchanalyzer); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
        r-dplyr=1.0.10: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
