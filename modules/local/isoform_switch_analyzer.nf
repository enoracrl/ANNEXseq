process ISOFORMSWITCHANALYZER {

    script:
    """
    IsoformSwitchAnalyzeR.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}
