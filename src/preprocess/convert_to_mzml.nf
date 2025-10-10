workflow convert_to_mzml {
    take:
    input_path

    main:
    if (params.is_timstof) {
        mzml_direct_convert = convert_bruker_d(input_path)
        mzml = adjust_mzML(mzml_direct_convert)
    } else {
        mzml = convert_thermo_raw(input_path)
    }
        
    emit:
    mzml
}

process convert_thermo_raw {
    cpus 2
    memory "8 GB"

    label 'msconvert_image'

	publishDir "${params.outdir}/mzmls", mode: 'copy', enabled: params.keep_mzmls

    input:
    path input_raw

    output:
    path "${input_raw.baseName}.mzML"

    script:
    """
    wine msconvert --mzML --zlib=off --filter "peakPicking vendor msLevel=1-" ${input_raw} --outfile ${input_raw.baseName}.mzML
    """
}

process convert_bruker_d {
    cpus { params.tdf2mzml_threads }
    memory "8 GB"
    
    label 'tdf2mzml_image'

    input:
    path input_d

    output:
    path "${input_d.baseName}.mzML"

    script:
    """
    export MKL_NUM_THREADS=${params.tdf2mzml_threads}
    export NUMEXPR_NUM_THREADS=${params.tdf2mzml_threads}
    export OMP_NUM_THREADS=${params.tdf2mzml_threads}

    tdf2mzml -i ${input_d} --compression "zlib" -o ${input_d.baseName}.mzML
    """
}

process adjust_mzML {
    cpus 2
    memory "8 GB"

    label 'msconvert_image'

	publishDir "${params.outdir}/mzmls", mode: 'copy', enabled: params.keep_mzmls

    input:
    path input_mzML

    output:
    path "uncompressed/${input_mzML.baseName}.mzML"

    script:
    """
    # msamanda needs explicit "scan=" in the id of a scan (not there in e.g. TimsTOF converted mzML data)
    sed -e 's/<spectrum\\(.*\\) id="index=\\([0-9]*\\)\\(.*\\)/<spectrum\\1 id="index=\\2 scan=\\2\\3/;s/spectrumRef="\\(.*\\)index=\\([0-9]*\\)\\(.*\\)"/spectrumRef="\\1index=\\2 scan=\\2\\3"/' ${input_mzML} > ${input_mzML.baseName}.reindex.mzML

    # after this, the mzML index needs to be removed, or just the compression???
    ## docker run --rm -it -v ./:/data/ proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25073-842baef wine msconvert --mzML -o ${input_mzML}.noindex_nocompress.mzML --noindex --zlib=off ${input_mzML};

    # uncompress the mzML with re-named id-tags (this re-calculates also the index, but for X!Tandem it is needed uncompressed anyways)
    mkdir uncompressed
    wine msconvert --mzML --zlib=off -o uncompressed --outfile ${input_mzML.baseName}.mzML ${input_mzML.baseName}.reindex.mzML
    
    # cleanup
    rm ${input_mzML.baseName}.reindex.mzML
    """
}


process split_mzml_into_chunks {
    cpus 2
    memory "8 GB"

    label 'msconvert_image'

    input:
    val chunksize
    path mzml_file

    output:
    tuple val(mzml_file.baseName), path("${mzml_file.baseName}--*.mzML")

    script:
    """
    CHUNKSIZE=${chunksize}

    # get only MS2 spectra
    wine msconvert --mzML --filter "msLevel 2-" "${mzml_file}" --outfile "${mzml_file.baseName}-MS2only.mzML"

    # get the last spectrum 
    MAX_INDEX=\$(tac ${mzml_file.baseName}-MS2only.mzML | grep -m1 "<spectrum.*index=" | sed 's;.*index="\\([0-9]*\\)".*;\\1;')
    echo "max index: \${MAX_INDEX}"

    for ((i=0; i<=MAX_INDEX; i+=CHUNKSIZE)); do
        START_INDEX=\$i
        END_INDEX=\$((i+CHUNKSIZE-1))
        if ((END_INDEX > MAX_INDEX)); then
            END_INDEX=\$MAX_INDEX
        fi
        echo "Processing spectra from index \${START_INDEX} to \${END_INDEX}"

        wine msconvert --mzML --filter "index [\${START_INDEX},\${END_INDEX}]" "${mzml_file.baseName}-MS2only.mzML" --outfile "${mzml_file.baseName}--\${START_INDEX}_\${END_INDEX}.mzML"
    done
    """
}
