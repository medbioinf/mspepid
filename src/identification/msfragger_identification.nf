include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

workflow msfragger_identification {
    take:
    default_params_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da
    // msfragger-specific runtime/configuration values passed from main.nf
    msfragger_threads
    msfragger_mem_gb
    msfragger_spectrum_id_pattern
    msfragger_scan_id_pattern
    msfragger_db_split
    msfragger_calibrate
    convert_psm_tsv_mem
    enhance_psm_tsv_mem
    use_only_rank1_psms
    ms2pip_model_dir
    ms2rescore_model
    ms2rescore_threads
    ms2rescore_mem
    ms2rescore_chunk_size
    oktoberfest_intensity_model
    oktoberfest_irt_model
    oktoberfest_memory
    oktoberfest_to_pin_memory
    oktoberfest_forks
    percolator_threads
    percolator_mem
    outdir

    main:
    fragger_params_file = adjust_msfragger_param_file(default_params_file, precursor_tol_ppm, fragment_tol_da, fasta)
    
    fragger_results = identification_with_msfragger(fasta, mzmls, fragger_params_file)
    fragger_results_pepxml = fragger_results.pepxml.flatten()
    
    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(fragger_results_pepxml, 'pepxml', 'msfragger', convert_psm_tsv_mem, enhance_psm_tsv_mem, outdir, use_only_rank1_psms)
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'msfragger', percolator_threads, percolator_mem, outdir)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.pepXML')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msfragger_spectrum_id_pattern, 'msfragger', ms2pip_model_dir, ms2rescore_model, ms2rescore_threads, ms2rescore_mem, ms2rescore_chunk_size, outdir)
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msfragger_scan_id_pattern, 'msfragger', fragment_tol_da, oktoberfest_intensity_model, oktoberfest_irt_model, oktoberfest_memory, oktoberfest_to_pin_memory, oktoberfest_forks, outdir)
    
    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'msfragger', percolator_threads, percolator_mem, outdir)
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'msfragger', percolator_threads, percolator_mem, outdir)
}


process adjust_msfragger_param_file {
    cpus 2
    memory "1 GB"

    label 'python_image'

    input:
    path fragger_params_file
    val precursor_tol_ppm
    val fragment_tol_da
    path fasta

    output:
    path "adjusted_fragger.params"

    script:
    """
    cp ${fragger_params_file} adjusted_fragger.params

    sed -i "s;^database_name.*;database_name = ${fasta};" adjusted_fragger.params
    sed -i "s;^num_threads.*;num_threads = ${msfragger_threads};" adjusted_fragger.params

    sed -i 's;precursor_mass_lower =.*;precursor_mass_lower = -${precursor_tol_ppm};' adjusted_fragger.params
    sed -i 's;precursor_mass_upper =.*;precursor_mass_upper = ${precursor_tol_ppm};' adjusted_fragger.params
    sed -i 's;precursor_mass_units =.*;precursor_mass_units = 1;' adjusted_fragger.params

    sed -i 's;fragment_mass_tolerance =.*;fragment_mass_tolerance = ${fragment_tol_da};' adjusted_fragger.params
    sed -i 's;fragment_mass_units =.*;fragment_mass_units = 0;' adjusted_fragger.params

    sed -i 's;calibrate_mass =.*;calibrate_mass = ${msfragger_calibrate};' adjusted_fragger.params

    sed -i "s;^decoy_prefix.*;decoy_prefix = DECOY_;" adjusted_fragger.params

    sed -i "s;^output_format.*;output_format = pepxml;" adjusted_fragger.params
    
    sed -i "s;^output_report_topN.*;output_report_topN = 5;" adjusted_fragger.params
    """
}


process identification_with_msfragger {
    cpus { msfragger_threads }
    memory { msfragger_mem_gb + " GB" }

    label 'msfragger_image'
    
    publishDir "${outdir}/msfragger", mode: 'copy'

    input:
    path fasta
    path mzmls
    path fragger_param_file

    output:
    path "*.pepXML", emit: pepxml

    script:
    if (msfragger_db_split > 0) {
        """
        msfragger_pep_split.py ${msfragger_db_split} "java -Xmx${msfragger_mem_gb}G -jar" /home/mambauser/MSFragger-4.2/MSFragger-4.2.jar ${fragger_param_file} ${mzmls} 
        """
    } else {
        """
        java -Xmx${msfragger_mem_gb}G -jar /home/mambauser/MSFragger-4.2/MSFragger-4.2.jar ${fragger_param_file} ${mzmls} 
        """
    }
}
