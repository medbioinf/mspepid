include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

// msamanda needs explicit "scan=" in the id of a scan (not there in e.g. TimsTOF converted mzML data)
// 1) ms-convert with "--noindex"
// 2) sed -i -e 's;<spectrum\(.*\) id="index=\([0-9]*\)\(.*\);<spectrum\1 id="index=\2 scan=\2\3;' filtered.mzML

// TODO:
// - params "LoadedProteinsAtOnce" adjustable
// - params "LoadedSpectraAtOnce" adjustable

/**
 * Executes the identification using MSAmanda
 */
workflow msamanda_identification {
    take:
    msamanda_config_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da
    // msamanda-specific runtime/config values passed from main.nf
    msamanda_threads
    msamanda_mem
    msamanda_spectrum_id_pattern
    msamanda_scan_id_pattern
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
    msamanda_results = identification_with_msamanda(msamanda_config_file, fasta, mzmls, precursor_tol_ppm, fragment_tol_da)

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(msamanda_results.msamanda_csv, 'msamanda', 'msamanda', convert_psm_tsv_mem, enhance_psm_tsv_mem, outdir, use_only_rank1_psms)
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'msamanda', percolator_threads, percolator_mem, outdir)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('_msamanda.csv')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msamanda_spectrum_id_pattern, 'msamanda', ms2pip_model_dir, ms2rescore_model, ms2rescore_threads, ms2rescore_mem, ms2rescore_chunk_size, outdir)
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msamanda_scan_id_pattern, 'msamanda', fragment_tol_da, oktoberfest_intensity_model, oktoberfest_irt_model, oktoberfest_memory, oktoberfest_to_pin_memory, oktoberfest_forks, outdir)

    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'msamanda', percolator_threads, percolator_mem, outdir)
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'msamanda', percolator_threads, percolator_mem, outdir)
}


process identification_with_msamanda {
    cpus { msamanda_threads }
    memory { msamanda_mem }

    label 'msamanda_image'

    publishDir "${outdir}/msamanda", mode: 'copy'

    input:
    path msamanda_config_file
    path fasta
    path mzmls
    val precursor_tol_ppm
    val fragment_tol_da

    output:
    path "${mzmls.baseName}_msamanda.csv", emit: msamanda_csv

    script:
    """
    cp ${msamanda_config_file} adjusted_msamanda_settings.xml
    sed -i 's;<MS1Tol[^<]*</MS1Tol>;<MS1Tol Unit="ppm">${precursor_tol_ppm}</MS1Tol>;' adjusted_msamanda_settings.xml
    sed -i 's;<MS2Tol[^<]*</MS2Tol>;<MS2Tol Unit="Da">${fragment_tol_da}</MS2Tol>;' adjusted_msamanda_settings.xml

    # MSAmanda command line arguments:
    # Required: -s spectrumFile     single .mgf or .mzml file, or folder with multiple .mgf and .mzml files
    # Required: -d proteinDatabase  single .fasta file or folder with multiple .fasta files, which will be combined into one
    # Required: -e settings.xml
    # Optional: -f fileformat       choose 1 for .csv and 2 for .mzid, default value is 1
    # Optional: -o outputfilename   file or folder where the output should be saved, default path is location of Spectrum file

    # create a short-name symlink or copy for the FASTA to work around MSAmanda's filename length limit (must be < 80 chars)
    if [ -e db.fasta ]; then rm -f db.fasta; fi
    ln -s "${fasta}" db.fasta || cp "${fasta}" db.fasta

    MSAmanda -s ${mzmls} -d db.fasta -e adjusted_msamanda_settings.xml -f 1 -o ${mzmls.baseName}_msamanda.csv
    """
}
