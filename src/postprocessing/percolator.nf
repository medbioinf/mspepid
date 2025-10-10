/**
 * Executes percolator for the given PIN files
 *
 * @return percolated PSMs (TSV pout file)
 */
workflow psm_percolator {
    take:
    pin_files
    searchengine

    main:
    pout_files = run_percolator(pin_files, searchengine)

    emit:
    pout_files
}


process run_percolator {
    cpus  { params.percolator_threads }
    memory { params.percolator_mem }

    label 'percolator_image'

	publishDir "${params.outdir}/${searchengine}", mode: 'copy'

    input:
    path pin_file
    val searchengine

    output:
    path "${pin_file.baseName}.pout"

    script:
    """
    percolator --num-threads ${params.percolator_threads} --only-psms --post-processing-tdc --search-input concatenated --results-psms ${pin_file.baseName}.pout ${pin_file} 
    """
}
