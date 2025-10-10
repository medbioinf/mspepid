workflow ms2rescore_workflow {
    take:
    psm_tsvs_and_mzmls
    psm_tsvs
    mzmls
    spectrum_id_pattern
    searchengine

    main:
    ms2pip_model_dir = Channel.fromPath(params.ms2pip_model_dir, type: 'dir').first()
    check_or_download_model(ms2pip_model_dir, params.ms2rescore_model)

    ms2rescore_pre_pins = run_chunked_ms2rescore(psm_tsvs_and_mzmls, psm_tsvs, mzmls, spectrum_id_pattern, params.fragment_tol_da, ms2pip_model_dir)
    ms2rescore_pins = correct_psm_utils_pins(ms2rescore_pre_pins, searchengine)

    emit:
    ms2rescore_pins
}


process run_chunked_ms2rescore {
    cpus  { params.ms2rescore_threads }
    memory { params.ms2rescore_mem }

    label 'python_image'

    input:
    tuple val(psm_utils_tsvs), val(mzml_for_psms)
    path psm_tsvs
    path mzmls
    val spectrum_id_pattern
    val fragment_tol_da
    path ms2pip_model_dir
    
    output:
    path "*.ms2rescore.pin", emit: features_file
    
    script:
    """
    chunked_ms2rescore.py -psms_file ${psm_utils_tsvs} \
        -spectra ${mzml_for_psms} \
        -model ${params.ms2rescore_model} -model_dir "${ms2pip_model_dir}" \
        -ms2_tolerance ${fragment_tol_da} \
        -spectrum_id_pattern '${spectrum_id_pattern}' \
        -processes ${params.ms2rescore_threads} \
        -chunk_size ${params.ms2rescore_chunk_size} \
        -out_file "${psm_utils_tsvs}.ms2rescore.pin"
    """
}


process correct_psm_utils_pins {
    cpus  2
    memory '8 GB'

    label 'python_image'

	publishDir "${params.outdir}/${searchengine}", mode: 'copy'

    input:
    path psm_utils_pins
    val searchengine

    output:
    path "${psm_utils_pins.baseName}.corrected.pin"

    script:
    """
    # correct the PIN file by moving the scan number to third column and adding correct SpecId (increasing integer)
    awk '{FS="\t";OFS="\t"; if (NR>1) { \$3=\$1; \$1=NR-1; gsub(".*=", "", \$3)  } print}' ${psm_utils_pins} > ${psm_utils_pins.baseName}.corrected.pin
    """
}


process check_or_download_model {
    cpus 1
    memory '2 GB'
    maxForks 1  // this makes sure that the download is only performed once, not more in parallel

    label 'python_image'

    input:
    path model_dir
    val ms2rescore_model

    script:
    """
    ms2rescore_check_or_download_model.py -ms2pip_model ${ms2rescore_model} -model_dir "${model_dir}"
    """
}