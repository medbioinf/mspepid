include {convert_and_enhance_psm_tsv} from '../postprocessing/convert_and_enhance_psm_tsv.nf'
include {psm_percolator; psm_percolator as ms2rescore_percolator; psm_percolator as oktoberfest_percolator} from '../postprocessing/percolator.nf'
include {ms2rescore_workflow} from '../postprocessing/ms2rescore.nf'
include {oktoberfest_rescore_workflow} from '../postprocessing/oktoberfest.nf'

/**
 * Executes the identification using Sage
 *
 * @return tuples containing the Sage results as [pin, tsv] files for each mzML
 */
workflow sage_identification {
    take:
    default_config_file
    fasta
    mzmls
    precursor_tol_ppm
    fragment_tol_da

    main:
    sage_config_file = adjust_sage_config(default_config_file, precursor_tol_ppm, fragment_tol_da)

    // all mzMLs are processed at once in sage for now -> much faster
    sage_results = identification_with_sage(sage_config_file, fasta, mzmls.collect())
    separated_results = separate_sage_results(sage_results.sage_tsv)

    psm_tsvs_and_pin = convert_and_enhance_psm_tsv(separated_results.sage_tsv.flatten(), 'sage_tsv', 'sage')
    psm_tsvs = psm_tsvs_and_pin.psm_tsv
    pin_files = psm_tsvs_and_pin.pin_file

    psm_percolator(pin_files, 'sage')

    psm_tsvs_and_mzmls = psm_tsvs.map { it -> [ it.name, it.name.take(it.name.lastIndexOf('.sage')) ] }
    ms2rescore_pins = ms2rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.sage_spectrum_id_pattern, 'sage')
    oktoberfest_pins = oktoberfest_rescore_workflow(psm_tsvs_and_mzmls, psm_tsvs.collect(), mzmls.collect(), params.sage_scan_id_pattern, 'sage')

    // perform percolation
    ms2rescore_percolator(ms2rescore_pins.ms2rescore_pins, 'sage')
    oktoberfest_percolator(oktoberfest_pins.oktoberfest_pins, 'sage')
}


process adjust_sage_config {
    cpus 2
    memory "1 GB"

    label 'python_image'

    input:
    path default_config_file
    val precursor_tol_ppm
    val fragment_tol_da

    output:
    path "adjusted_sage_config.json"

    script:
    """
#!/usr/bin/env python
import json

# Opening JSON file
with open("${default_config_file}", 'r') as openfile:
    # Reading from json file
    json_object = json.load(openfile)

# adjust the tolerances
json_object["precursor_tol"] = {'ppm': [-${precursor_tol_ppm}, ${precursor_tol_ppm}]}
json_object["fragment_tol"] = {'da': [-${fragment_tol_da}, ${fragment_tol_da}]}

# adjust database prefilter
json_object["database"]["prefilter_chunk_size"] = ${params.sage_prefilter_chunk_size}
json_object["database"]["prefilter"] = (str("${params.sage_prefilter}").strip().lower() == "true")
json_object["database"]["prefilter_low_memory"] = False

# Writing to sample.json
with open("./adjusted_sage_config.json", "w") as outfile:
    json.dump(json_object, outfile)
    """
}


process identification_with_sage {
    cpus { params.sage_threads }
    memory { params.sage_mem }

    label 'sage_image'

    input:
    path sage_config_file
    path fasta
    path mzmls

    output:
    path "results.sage.tsv", emit: sage_tsv

    script:
    """
    RAYON_NUM_THREADS=${params.sage_threads} sage ${sage_config_file} -f ${fasta} --batch-size ${params.sage_threads} --write-pin --output_directory ./ ${mzmls}
    """
}


process separate_sage_results {
    cpus 2
    memory "1 GB"

    label 'python_image'

    publishDir "${params.outdir}/sage", mode: 'copy'

    input:
    path sage_tsv

    output:
    path "*.sage.tsv", emit: sage_tsv

    script:
    """
    # process the tsv file and create one file for each input file
    for filename in \$(awk 'NR>1{a[\$5]++} END{for(b in a) print b}' ${sage_tsv});
    do
        head -n1 ${sage_tsv} > \${filename}.sage.tsv
        awk -v f1="\${filename}" '\$5==f1' ${sage_tsv} >> \${filename}.sage.tsv
    done
    """
}
