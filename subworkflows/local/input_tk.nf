/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_TK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    input_path

    main:
    /*
     * Check samplesheet is valid
     */
    SAMPLESHEET_CHECK ( samplesheet, input_path )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it, params.genomes) }
        .map { it -> [ it[0], it[8] ] } // a modif
        .set { ch_model }

    emit:
    ch_model // [ sample, model]
}


def get_sample_info(LinkedHashMap sample, LinkedHashMap genomeMap) {
    def meta = [:]
    meta.id  = sample.sample

    // Resolve fasta and gtf file if using iGenomes
    def fasta = false
    def gtf   = false
    if (sample.fasta) {
        if (genomeMap && genomeMap.containsKey(sample.fasta)) {
            fasta = file(genomeMap[sample.fasta].fasta, checkIfExists: true)
            gtf   = file(genomeMap[sample.fasta].gtf, checkIfExists: true)
        } else {
            fasta = file(sample.fasta, checkIfExists: true)
        }
    }

    // Check if input file and gtf file exists
    input_file = sample.input_file ? file(sample.input_file, checkIfExists: true) : null
    gtf        = sample.gtf        ? file(sample.gtf, checkIfExists: true)        : gtf

    return [ meta, input_file, sample.barcode, fasta, gtf, sample.is_transcripts.toBoolean(), fasta.toString()+';'+gtf.toString(), sample.nanopolish_fast5 ]
}
