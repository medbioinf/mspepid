/**
 * Adds decoys and/or entapments to the FASTA file.
 *
 * @param fasta Path to the FASTA file
 *
 * @return Path to the new FASTA file
 */
workflow create_entrapment_database {
    take:
        fasta
        fold
        fdrbench_mem_gb

    main:
        entrapment_fasta = call_entrapment_database(fasta, fold, fdrbench_mem_gb)
        
    emit:
        entrapment_fasta
}


/**
 * Adds entrapments to the FASTA file using FDRBench.
 * https://doi.org/10.1101/2024.06.01.596967
 *
 * @param fasta Path to the FASTA file
 * @param fold Fold change for entrapment
 *
 * @return Path to the new FASTA file
 */
process call_entrapment_database {
    cpus 1
    memory "${ memory_limit }.GB" //fdrbench_mem_gb + " GB"

    label 'fdrbench_image'

    input: 
    path fasta
    val fold
    val memory_limit   //fdrbench_mem_gb + " GB"

    output:
    path "${fasta.baseName}-entrapment.fasta"

    script:
    """
    java -Xmx${memory_limit}G -jar /opt/fdrbench/fdrbench.jar -db ${fasta} -o ${fasta.baseName}-entrapment.fasta -fold ${fold} -level protein -entrapment_label ENTRAPMENT_ -entrapment_pos 0 -uniprot -check
    # 'Reheader' to add entrapment index to database and accession part of the header
    # and remove empty entrapment sequences (which can appear if the original sequence has many Xs)
    # The following sed command performs two operations:
    # 1. Substitutes FASTA headers to include the entrapment index in both the database and accession parts.
    #    Regex breakdown:
    #      ^>ENTRAPMENT_(.+)\\|(.+)\\|(.+)_([0-9]+)\$
    #        - Matches headers starting with '>ENTRAPMENT_' followed by three fields separated by '|', with the last field ending in '_[number]'.
    #      >ENTRAPMENT_\\4_\\1|ENTRAPMENT_\\4_\\2|\\3_\\4
    #        - Rewrites the header to include the entrapment index (\\4) in both the database and accession parts.
    # 2. Removes empty entrapment sequences (headers followed by an empty line).
    #    Control flow:
    #      \$!N;/>.*\\n\$/d;P;D
    #        - Reads two lines at a time; if a header is followed by an empty line, deletes both.

    sed -r -i -e "s;^>ENTRAPMENT_(.+)\\|(.+)\\|(.+)_([0-9]+)\$;>ENTRAPMENT_\\4_\\1|ENTRAPMENT_\\4_\\2|\\3_\\4;g" -e '\$!N;/>.*\\n\$/d;P;D'  ${fasta.baseName}-entrapment.fasta
    """
}
