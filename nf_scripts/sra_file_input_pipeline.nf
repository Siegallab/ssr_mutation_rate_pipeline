/*
 * Workflow for running variant calling starting with a list of SRA ids
 */

nextflow.enable.dsl=2

include { align_call_join_recall } from './generic_workflow.nf'

process get_SRA {

    // uses sra_tools to get files for SRA_id
    // if sra_tools prefetch fails (happens seemingly randomly for some 
    // files, tries to download the files manually, assuming they are
    // .fastq.gz)

    publishDir "${params.out_dir}/${params.sra_out}", mode:'copy', overwrite: true

    input: val(SRA_id)

    output: tuple val(SRA_id), file("${SRA_id}_pass*.fastq.gz")

    script:
    """
    module load ${params.SRA_TOOLS}
    # note that prefetch randomly fails for some files
    # if this happens, get link and download manually
    if prefetch $SRA_id 2> prefetch_error_output.txt; then
        fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ${SRA_id}/${SRA_id}.sra
    else
        # remove directory that just got created
        rm -r $SRA_id
        # get link to folder containing files prefetch was trying to get
        #folder_link=`tail -n 1 prefetch_error_output.txt | sed -n -e 's/^.*https:/https:/p' | xargs dirname`
        folder_link=`tail -n 1 prefetch_error_output.txt | sed -n -e 's/^.*https:/https:/p'`
        wget --recursive --no-parent -nH --cut-dirs=4 \$folder_link
        fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip ${SRA_id}
    fi
    """
}

workflow {
	Channel.fromList(params.id_dict) \
    	| map{it[0]} \
    	| get_SRA \
		| align_call_join_recall
}