include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

/**
 * Exports the identification using Comet configured by a SDRF files
 */
workflow comet_identification {
    take:
    default_params_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da

    main:
    comet_params_file = adjust_comet_param_file(default_params_file, precursor_tol_ppm, fragment_tol_da)
    
    comet_mzids = identification_with_comet(fasta, mzmls, comet_params_file)
    comet_mzids = comet_mzids.flatten()
    
    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(comet_mzids, 'mzid', 'comet')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'comet')

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.comet_spectrum_id_pattern, 'comet')
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.comet_scan_id_pattern, 'comet')
    
    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'comet')
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'comet')
}


process adjust_comet_param_file {
    cpus 2
    memory "1 GB"

    label 'python_image'

    input:
    path comet_params_file
    val precursor_tol_ppm
    val fragment_tol_da

    output:
    path "adjusted_comet.params"

    script:
    """
    cp ${comet_params_file} adjusted_comet.params
    
    sed -i 's;peptide_mass_tolerance_upper =.*;peptide_mass_tolerance_upper = ${precursor_tol_ppm};' adjusted_comet.params
    sed -i 's;peptide_mass_tolerance_lower =.*;peptide_mass_tolerance_lower = -${precursor_tol_ppm};' adjusted_comet.params
    sed -i 's;peptide_mass_units =.*;peptide_mass_units = 2;' adjusted_comet.params

    sed -i 's;fragment_bin_tol =.*;fragment_bin_tol = ${fragment_tol_da};' adjusted_comet.params
    
    sed -i "s;^num_threads.*;num_threads = ${params.comet_threads};" adjusted_comet.params

    sed -i "s;^output_sqtfile.*;output_sqtfile = 0;" adjusted_comet.params
    sed -i "s;^output_txtfile.*;output_txtfile = 0;" adjusted_comet.params
    sed -i "s;^output_pepxmlfile.*;output_pepxmlfile = 0;" adjusted_comet.params
    sed -i "s;^output_mzidentmlfile.*;output_mzidentmlfile = 1;" adjusted_comet.params
    sed -i "s;^output_percolatorfile.*;output_percolatorfile = 0;" adjusted_comet.params
    
    sed -i "s;^num_output_lines.*;num_output_lines = 5;" adjusted_comet.params
    """
}


process identification_with_comet {
    cpus { params.comet_threads }
    memory { params.comet_mem }

    label 'comet_image'

	publishDir "${params.outdir}/comet", mode: 'copy'

    input:
    path fasta
    path mzmls
    path comet_param_file

    output:
    path "*.mzid"

    script:
    """
    comet -P${comet_param_file} -D${fasta} ${mzmls}
    """
}
