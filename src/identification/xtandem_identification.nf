include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

/**
 * Exports the identification using Comet configured by a SDRF files
 */
workflow xtandem_identification {
    take:
    xtandem_config_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da
    // xtandem-specific runtime/configuration values passed from main.nf
    xtandem_spectrum_id_pattern
    xtandem_scan_id_pattern
    xtandem_threads
    xtandem_mem
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
    (xtandem_param_files, taxonomy_file) = create_xtandem_params_files_from_default(xtandem_config_file, fasta, mzmls, precursor_tol_ppm, fragment_tol_da)
    xtandem_param_files = xtandem_param_files.flatten()

    tandem_xmls = identification_with_xtandem(xtandem_param_files, taxonomy_file, fasta, mzmls.collect())
    tandem_xmls = tandem_xmls.flatten()

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(tandem_xmls, 'xtandem', 'xtandem', convert_psm_tsv_mem, enhance_psm_tsv_mem, outdir, use_only_rank1_psms)
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'xtandem', percolator_threads, percolator_mem, outdir)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.xtandem_identification')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), xtandem_spectrum_id_pattern, 'xtandem', ms2pip_model_dir, ms2rescore_model, ms2rescore_threads, ms2rescore_mem, ms2rescore_chunk_size, outdir)
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), xtandem_scan_id_pattern, 'xtandem', fragment_tol_da, oktoberfest_intensity_model, oktoberfest_irt_model, oktoberfest_memory, oktoberfest_to_pin_memory, oktoberfest_forks, outdir)

    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'xtandem', percolator_threads, percolator_mem, outdir)
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'xtandem', percolator_threads, percolator_mem, outdir)
}

/**
 * Creates a X!Tandem params file from the given default file for the mzML files
 * @param xtandem_config_file The default config file
 * @param sdrf The FASTA file
 * @param max_missed_clavages maximum number of missed cleavages
 * @param max_parent_charge  maximum parent charge

 * @return The XTandem params for each file in the SDRF and the according taxonomy file
 */
process create_xtandem_params_files_from_default {
    cpus 2
    memory "1 GB"

    label 'python_image'

    input:
    path xtandem_config_file
    path fasta
    path mzmls
    val precursor_tol_ppm
    val fragment_tol_da

    output:
    path "*.xtandem_input.xml"
    path "xtandem_taxonomy.xml"

    script:
    """
    # write the taxonomy file
    echo '<?xml version="1.0"?>
<bioml label="x! taxon-to-file matching list">
  <taxon label="sample_species">
    <file format="peptide" URL="${fasta}" />
  </taxon>
</bioml>' > xtandem_taxonomy.xml

    # adjust parameters in the default file
    cp ${xtandem_config_file} ${mzmls.baseName}.xtandem_input.xml

    sed -i 's;<note type="input" label="list path, taxonomy information">[^<]*</note>;<note type="input" label="list path, taxonomy information">xtandem_taxonomy.xml</note>;' ${mzmls.baseName}.xtandem_input.xml

    sed -i 's;<note type="input" label="spectrum, path">[^<]*</note>;<note type="input" label="spectrum, path">${mzmls}</note>;' ${mzmls.baseName}.xtandem_input.xml
    sed -i 's;<note type="input" label="output, path">[^<]*</note>;<note type="input" label="output, path">${mzmls.baseName}.xtandem_identification.t.xml</note>;' ${mzmls.baseName}.xtandem_input.xml

    sed -i 's;<note type="input" label="spectrum, fragment monoisotopic mass error">[^<]*</note>;<note type="input" label="spectrum, fragment monoisotopic mass error">${fragment_tol_da}</note>;' ${mzmls.baseName}.xtandem_input.xml
    sed -i 's;<note type="input" label="spectrum, fragment monoisotopic mass error units">[^<]*</note>;<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>;' ${mzmls.baseName}.xtandem_input.xml

    sed -i 's;<note type="input" label="spectrum, parent monoisotopic mass error minus">[^<]*</note>;<note type="input" label="spectrum, parent monoisotopic mass error minus">${precursor_tol_ppm}</note>;' ${mzmls.baseName}.xtandem_input.xml
    sed -i 's;<note type="input" label="spectrum, parent monoisotopic mass error plus">[^<]*</note>;<note type="input" label="spectrum, parent monoisotopic mass error plus">${precursor_tol_ppm}</note>;' ${mzmls.baseName}.xtandem_input.xml
    sed -i 's;<note type="input" label="spectrum, parent monoisotopic mass error units">[^<]*</note>;<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>;' ${mzmls.baseName}.xtandem_input.xml

    sed -i 's;<note type="input" label="spectrum, threads">[^<]*</note>;<note type="input" label="spectrum, threads">${xtandem_threads}</note>;' ${mzmls.baseName}.xtandem_input.xml

    # rename absolute paths to current path, to allow for clean passing on in workflow
    workDir=\$(pwd)
    for file in *.xtandem_input.xml; do
        sed -i "s;\$workDir/;;g" \$file
    done

    #sed -i "s;\$workDir/;;g" xtandem_taxonomy.xml
    """
}

/**
 * Performs the identifications with XTandem
 */
process identification_with_xtandem {
    cpus { xtandem_threads }
    memory { xtandem_mem }

    label 'xtandem_image'
    
    publishDir "${outdir}/xtandem", mode: 'copy'

    input:
    path xtandem_param_file
    path taxonomy_file
    path fasta
    path mzmls

    output:
    path "*.t.xml"

    script:
    """
    tandem $xtandem_param_file
    """
}
