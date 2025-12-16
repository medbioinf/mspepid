/**
 * Runs oktoberfest rescoring for the given PSMs and mzML files.
 * 
 * @param psm_tsvs_and_mzmls: A tuple containing the PSM utils TSV files and the mzML files for the PSMs.
 * @param psm_tsvs: The PSM TSV files.
 * @param mzmls: The mzML files. 
 * @param scan_id_regex: A regex pattern to extract the scan number from the spectrum ID.
 *
 * @return: The oktoberfest rescored PSMs in TSV format.
 */
workflow oktoberfest_rescore_workflow {
    take:
    psm_tsvs_and_mzmls
    psm_tsvs
    mzmls
    scan_id_regex
    searchengine
    fragment_tol_da
    oktoberfest_intensity_model
    oktoberfest_irt_model
    oktoberfest_memory
    oktoberfest_to_pin_memory
    oktoberfest_forks
    outdir

    main:
    oktoberfest_features = run_oktoberfest_feature_gen(psm_tsvs_and_mzmls, psm_tsvs, mzmls, fragment_tol_da, scan_id_regex, oktoberfest_intensity_model, oktoberfest_irt_model, oktoberfest_memory, oktoberfest_forks)
    oktoberfest_pins = oktoberfest_features_to_pin(oktoberfest_features, searchengine, oktoberfest_to_pin_memory, outdir)


    emit:
    oktoberfest_pins
}

/**
 * @param psm_tsvs_and_mzmls: A tuple containing the PSM utils TSV files and the mzML files for the PSMs.
 * @param psm_tsvs: The PSM TSV files.
 * @param mzmls: The mzML files. 
 * @param fragment_tol_da: The fragment tolerance for the rescoring.
 * @param scan_id_regex: A regex pattern to extract the scan number from the spectrum ID.
 * 
 * @return The oktoberfest rescored PSMs in TSV format.
 */
process run_oktoberfest_feature_gen {
    cpus 1
    maxForks { oktoberfest_forks }
    memory { oktoberfest_memory }

    label 'oktoberfest_image'

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val fragment_tol_da
    val scan_id_regex
    val oktoberfest_intensity_model
    val oktoberfest_irt_model
    
    
    output:
    path "${psm_utils_tsvs}.features.tsv"
    
    script:
    """
    oktoberfest_feature_gen.py \
        -out-folder ./oktoberfest_out \
        -psms-file ${psm_utils_tsvs} \
        -spectra-file ${mzml_for_psms} \
        -intensity-model ${oktoberfest_intensity_model} \
        -irt-model ${oktoberfest_irt_model} \
        -mass-tolerance ${fragment_tol_da} \
        -mass-tolerance-unit da \
        -scan-id-regex '${scan_id_regex}' \

    mv ./oktoberfest_out/results/none/rescore.tab "${psm_utils_tsvs}.features.tsv"

    # Clean up the output directory
    rm -r oktoberfest_out
    """
}

/**
 * @param okt_features_tsv: Oktoberfest feature file.
 * 
 * @return Oktoberfest feature file in PIN format ready to use with percolator.
 */
process oktoberfest_features_to_pin {
    cpus 1
    memory { oktoberfest_to_pin_memory }

    label 'oktoberfest_image'

	publishDir "${outdir}/${searchengine}", mode: 'copy'

    input:
    path okt_features_tsv
    val searchengine
    val oktoberfest_to_pin_memory
    val outdir

    output:
    path "${okt_features_tsv.baseName}.oktoberfest.pin"

    script:
    """
    oktoberfest_feature_to_pin.py \
        -in-file ${okt_features_tsv} \
        -out-file ./${okt_features_tsv.baseName}.oktoberfest.pin
    """
}
