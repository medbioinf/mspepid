include {convert_chunked_result_to_psm_utils; enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

include {split_mzml_into_chunks} from '../preprocess/convert_to_mzml.nf'

/**
 * Executes the identification using MS-GF+
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow msgfplus_identification {
    take:
    msgfplus_params_file
    fasta
    mzmls
    precursor_tol_ppm
    // msgfplus-specific runtime/configuration values passed from main.nf
    msgfplus_split_fasta
    msgfplus_split_input
    msgfplus_spectrum_id_pattern
    msgfplus_scan_id_pattern
    msgfplus_threads
    msgfplus_mem_gb
    msgfplus_instrument
    msgfplus_tasks
    msgfplus_merge_mem_gb
    msgfplus_mzid_mem_gb
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
    if (msgfplus_split_fasta > 0) {
        fasta = split_fasta(fasta)
        fasta = fasta.flatten()
    }
    fasta_index = build_msgfplus_index(fasta)

    if (msgfplus_split_input > 0) {
        chunked_mzmls = split_mzml_into_chunks(msgfplus_split_input, mzmls)
        mzmls_to_chunks = chunked_mzmls.transpose()
    } else {
        mzmls_to_chunks = mzmls.map{ it -> [it.baseName, it] }
    }

    fasta_idx_mzml_chunk_combo = fasta_index.combine(mzmls_to_chunks)

    // publish_results is true when we did not split fasta (single output per mzML)
    publish_results = (msgfplus_split_fasta == 0)
    msgfplus_results = identification_with_msgfplus(msgfplus_params_file, fasta_idx_mzml_chunk_combo, precursor_tol_ppm, publish_results)

    if ((msgfplus_split_fasta > 0)) {
        grouped_results = msgfplus_results.map { it ->
            tuple(
                it[0],  // original mzml basename
                it[1].name.take(it[1].name.lastIndexOf('-split')),  // mzML split
                it[1].name.substring(it[1].name.lastIndexOf('-split')+7, it[1].name.lastIndexOf('.mzid')), // FASTA split number
                it[1]  // mzid file
            )
        }.groupTuple(by: 1).map { it ->
            tuple(
                it[0][0], // original mzml basename
                it[1],    // mzML split
                it[3]     // collect all mzid files for this group
            )
        }
        fasta_merged_results = mzid_merger(grouped_results)
    } else {
        fasta_merged_results = msgfplus_results
    }

    psm_tsvs_with_mzml = convert_chunked_result_to_psm_utils(fasta_merged_results, 'mzid')

    grouped_results = psm_tsvs_with_mzml.groupTuple(by: 0)
    if (msgfplus_split_input > 0) {
        merged_results = merge_psms(grouped_results)
    } else {
        merged_results = psm_tsvs_with_mzml.map{ it -> it[1] }
    }

    psm_tsvs_and_pin = enhance_psm_tsv(merged_results, 'msgfplus', convert_psm_tsv_mem, enhance_psm_tsv_mem, outdir, use_only_rank1_psms)

    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'msgfplus', percolator_threads, percolator_mem, outdir)

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.mzid')) + '.mzML'  ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msgfplus_spectrum_id_pattern, 'msgfplus', ms2pip_model_dir, ms2rescore_model, ms2rescore_threads, ms2rescore_mem, ms2rescore_chunk_size, outdir)
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), msgfplus_scan_id_pattern, 'msgfplus', fragment_tol_da, oktoberfest_intensity_model, oktoberfest_irt_model, oktoberfest_memory, oktoberfest_to_pin_memory, oktoberfest_forks, outdir)

    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'msgfplus', percolator_threads, percolator_mem, outdir)
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'msgfplus', percolator_threads, percolator_mem, outdir)
}

process identification_with_msgfplus {
    cpus { msgfplus_threads }
    memory { msgfplus_mem_gb + " GB" }

    label 'msgfplus_image'

    publishDir "${outdir}/msgfplus", mode: 'copy', enabled: { publish_results }

    input:
    path msgfplus_params_file
    tuple path(fasta), path(canno), path(cnlcp), path(csarr), path(cseq), val(original_mzml_basename), path(mzml)
    val precursor_tol_ppm
    val publish_results

    output:
    tuple val(original_mzml_basename), path("${mzml.baseName}*.mzid")

    script:
    """
    cp ${msgfplus_params_file} adjusted_MSGFPlus_Params.txt
    sed -i 's;^PrecursorMassTolerance=.*;PrecursorMassTolerance=${precursor_tol_ppm};' adjusted_MSGFPlus_Params.txt
    sed -i 's;^InstrumentID=.*;InstrumentID=${msgfplus_instrument};' adjusted_MSGFPlus_Params.txt

    java -Xmx${msgfplus_mem_gb}G -jar /opt/msgfplus/MSGFPlus.jar -conf adjusted_MSGFPlus_Params.txt -s ${mzml} -d ${fasta} -thread ${msgfplus_threads} -tasks ${msgfplus_tasks} -o ${mzml.baseName}.mzid

    if [[ ${fasta} == *"-split"* ]]; then
        splitnum=\$(echo "${fasta}" | sed "s;.*-split-\\([0-9]*\\).fasta;\\1;")
        echo "renaming ${mzml.baseName}.mzid to ${mzml.baseName}-split-\${splitnum}.mzid"
        mv ${mzml.baseName}.mzid ${mzml.baseName}-split-\${splitnum}.mzid
    fi
    """
}


process split_fasta {
    cpus 2
    memory "8 GB"

    label 'python_image'

    input:
    path fasta

    output:
    path "${fasta.baseName}-split*.fasta", emit: fasta_parts

    script:
    """
    split_fasta.py -in_file ${fasta} -out_file_base ${fasta.baseName}-split -splits ${msgfplus_split_fasta} 
    """
}


process build_msgfplus_index {
    cpus { msgfplus_threads }
    memory { msgfplus_mem_gb + " GB" }

    label 'msgfplus_image'

    input:
    path fasta

    output:
    tuple path("${fasta}"), path("${fasta.baseName}.canno"), path("${fasta.baseName}.cnlcp"), path("${fasta.baseName}.csarr"), path("${fasta.baseName}.cseq")

    script:
    """
    java -Xmx${msgfplus_mem_gb}G -cp /opt/msgfplus/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d ${fasta} -tda 0 -o ./ -decoy DECOY_
    """
}


process merge_psms {
    cpus 2
    memory { msgfplus_merge_mem_gb + " GB" }

    label 'python_image'

    input:
    tuple val(original_mzml_basename), path(psm_tsvs)

    output:
    path "${original_mzml_basename}.mzid.psm_utils.tsv"

    script:
    """
    merge_chunked_psm_files.py --org_filebase ${original_mzml_basename} --out_filename ${original_mzml_basename}.mzid.psm_utils.tsv --files ${psm_tsvs}
    """
}


process mzid_merger {
    cpus 2
    memory { msgfplus_mzid_mem_gb + " GB" }

    label 'mzidmerger_image'

    publishDir "${outdir}/msgfplus", mode: 'copy'

    input:
    tuple val(original_mzml_basename), val(mzml_split), path(mzid_files)

    output:
    tuple val(original_mzml_basename), path("${mzml_split}.mzid")

    script:
    """
    mono /home/mambauser/mzidmerger/MzidMerger.exe -InDir "./" -Filter "*.mzid" -Out ${mzml_split}.merged.mzid
    mv ${mzml_split}.merged.mzid ${mzml_split}.mzid
    """
}