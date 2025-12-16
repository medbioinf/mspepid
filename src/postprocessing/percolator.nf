/**
 * Executes percolator for the given PIN files
 *
 * @return percolated PSMs (TSV pout file)
 */
workflow psm_percolator {
    take:
    pin_files
    searchengine
    percolator_threads
    percolator_mem
    outdir

    main:
    pout_files = run_percolator(pin_files, searchengine, percolator_threads, percolator_mem, outdir)

    emit:
    pout_files
}


process run_percolator {
    cpus  { percolator_threads }
    memory { percolator_mem }

    label 'percolator_image'

    publishDir "${outdir}/${searchengine}", mode: 'copy'

    input:
    path pin_file
    val searchengine
    val percolator_threads
    val percolator_mem
    val outdir

    output:
    path "${pin_file.baseName}.pout"

    script:
    """
    percolator --num-threads ${percolator_threads} --only-psms --post-processing-tdc --search-input concatenated --results-psms ${pin_file.baseName}.pout ${pin_file}
    """
}
