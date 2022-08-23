/*
 *
 * Manually remove duplicate reads and set up combined reference
 *
 * Manually removing duplicates is important because we merge two 
 * fasta files, so we want to avoid identical reads from 
 * the different files counting as 'duplicates'; to avoid this,
 * we will later run FreeBayes and HipSTR without ignoring duplicates
 *
 */

process subset_fastq_pe {

    // subset fastq files

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path(read_list)

    script:
    reads_in_1 = reads[0]
    reads_in_2 = reads[1]
    reads_out_1 = "subset_${pair_id}_${reads_in_1}"
    reads_out_2 = "subset_${pair_id}_${reads_in_2}"
    read_list = [reads_out_1, reads_out_2]
    """
    # use same seed for subsetting all sequences
    module load ${params.SAMTOOLS}
    ${params.seqtk} sample -s 123 ${reads[0]} ${params.subset_read_num} | bgzip > $reads_out_1
    ${params.seqtk} sample -s 123 ${reads[1]} ${params.subset_read_num} | bgzip > $reads_out_2
    """

}

process subset_fastq_se {

    // subset fastq files

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path(reads_out)

    script:
    reads_out = "subset_${pair_id}_${reads}"
    """
    # use same seed for subsetting all sequences
    ${params.seqtk} sample -s 123 ${reads} ${params.subset_read_num} > $reads_out
    """

}

process trim_pe {

    // perform trimming on paired-end reads

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("trimmed_{1,2}.fq.gz")

    script:
    """
    module load ${params.TRIMMOMATIC}
    java -jar \$TRIMMOMATIC_JAR \
    PE \
    -phred33 \
    -threads ${task.cpus} \
    ${reads[0]} \
    ${reads[1]} \
    trimmed_1.fq.gz \
    unpair_trimmed_1.fq.gz \
    trimmed_2.fq.gz \
    unpair_trimmed_2.fq.gz \
    ILLUMINACLIP:\$TRIMMOMATIC_HOME/{params.adapters}:2:30:10:8:true \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process trim_se {

    // perform trimming on single-end reads

    input:
    tuple val(id), path(read)

    output:
    tuple val(id), path("trimmed.fq.gz")

    script:
    """
    module load ${params.TRIMMOMATIC}
    java -jar \$TRIMMOMATIC_JAR \
    SE \
    -phred33 \
    -threads ${task.cpus} \
    $read \
    trimmed.fq.gz \
    ILLUMINACLIP:\$TRIMMOMATIC_HOME/{params.adapters}:2:30:10:8:true \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process idx_for_align {

    // set up indices of reference fasta for alignment

    input:
    path(ref)

    output:
    tuple path(ref), path("${ref_basename}.*")

    script:
    ref_basename = ref.getBaseName()
    """
    module load ${params.BWA}
    # reindex to have index files in current dir (fast)
    bwa index $ref
    """

}

process align_mem {

    input:
    tuple path(ref), path(ref_idx_files), val(id), path(reads)
     
    output:
    tuple val(id), path(out_sam)
    
    script:
    out_sam = "${id}_mem_aligned_reads.sam"
    read = reads.join(' ')
    readGroup = "@RG\\tID:${id}\\tLB:${id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${id}"
    """
    module load ${params.BWA}
    bwa mem \
    -O 5 \
    -w 500 \
    -K 100000000 \
    -v 3 \
    -t ${task.cpus} \
    -Y \
    -R \"${readGroup}\" \
    $ref \
    $read \
    > $out_sam
    """
}
        
process remove_dups {

    input:
    tuple val(id), path(aligned_reads)

    output:
    tuple val(id), path(out_dedup)

    script:
    out_dedup = "${id}_no_dups.bam"
    """
    module load ${params.GATK}
    gatk MarkDuplicatesSpark \
        -I $aligned_reads \
        -O $out_dedup \
        --remove-all-duplicates
    """

}

process merge_dedup_bams {

    input:
    tuple val(id), path(merge_input_file_list)

    output:
    tuple val(id), path(merge_output)

    script:
    merge_input = merge_input_file_list.join(' I=')
    prelim_merge_output = "${id}_no_dups_wrong_rg.bam"
    merge_output = "${id}_no_dups_merged.bam"
    """
    module load ${params.PICARD}
    java -jar \$PICARD_JAR \
        MergeSamFiles \
        I=${merge_input} \
        O=${prelim_merge_output} \
        VALIDATION_STRINGENCY=LENIENT
    # combine read groups
    java -jar \$PICARD_JAR \
        AddOrReplaceReadGroups \
        I=${prelim_merge_output} \
        O=$merge_output \
        RGID=${id} \
        RGLB=${id} \
        RGPL=${params.pl} \
        RGPU=${id} \
        RGPM=${params.pm} \
        RGSM=${id} \
        VALIDATION_STRINGENCY=LENIENT
    """

}

process sort {

    // Sort sam/bam files

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val publish_dir
    tuple val(id), path(aligned_reads)

    output:
    tuple val(id), path(out_bam)

    script:
    out_bam = "${id}_sorted.bam"
    prelim_out_bam = "${id}_prelim_sorted.bam"
    """
    module load ${params.PICARD}
    module load ${params.SAMTOOLS}
    # there's a PICARD bug that causes file to be misformed sometimes
    # get around this by looping until non-zero filesize
    # only works if input is a bam
    target_file_size=`samtools view $aligned_reads | wc -l`
    file_size=0
    while [ \$file_size -lt \$target_file_size ]
    do
        java -jar \$PICARD_JAR \
            SortSam \
            I=$aligned_reads \
            O=$prelim_out_bam \
            SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=LENIENT
        samtools quickcheck -q $prelim_out_bam
        quickcheck_err_stat=\$?
        if [ \$quickcheck_err_stat -ne 0 ]
        then
            file_size=0
        else
            file_size=`samtools view $prelim_out_bam | wc -l`
        fi
    done
    mv $prelim_out_bam $out_bam
    """ 
}

process get_low_coverage {

    // output .bed file for every sample containing regions in which coverage is below params.min_coverage

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val(publish_dir)
    tuple val(id), path(sorted_bam_file)

    output:
    path(coverage_file)

    script:
    coverage_file = "${id}_coverage.bed"
    """
    module load ${params.BEDTOOLS}
    bedtools genomecov -ibam $sorted_bam_file -bga -max ${params.min_coverage} | \
        awk '\$4 < ${params.min_coverage}' | \
        bedtools merge -i stdin -d 1 > \
        ${coverage_file}
    """
}

/*
 *
 * Read pairing-agnostic workflows that rely on a supplied params.PE bool
 *
 */

workflow subset_fastq {
    take:
        input_file_ch // [id, read(s)]

    main:
        if( params.PE )
            subset_fastq_pe(input_file_ch) | set{ subset_out_ch }
        else
            subset_fastq_se(input_file_ch) | set{ subset_out_ch }

    emit:
        subset_out_ch
}

workflow trim {
    take:
        input_file_ch // [id, read(s)]

    main:
        if( params.PE )
            trim_pe(input_file_ch) | set{ trim_out_ch }
        else
            trim_se(input_file_ch) | set{ trim_out_ch }

    emit:
        trim_out_ch
}

workflow perform_dedup {

    // perform deduplication

    take:
        align_out // [id, bam/sam]

    main:
        blank_ch = channel.value("")
        sort(blank_ch, align_out) | remove_dups

    emit:
        remove_dups.out

}

workflow postalign_process {

    // perform post-alignment duplicate removal and merging based on id

    take:
        publishdir_ch
            // Value channel containing directory for alignments;
            // if empty, alignments not saved
        align_out // [id, bam]

    main:
        dedup_ch = perform_dedup(align_out)
        // replace IDs based on id_dict, and merge data from any 
        // strains that need merging
        // This tries to merge even single bams, renames the ID in the 
        // read groups in the process
        id_dict_ch = channel.fromList(params.id_dict)
        id_dict_ch \
                // [ori_id, new_id]
            | join(dedup_ch, by: 0) \
                // [ori_id, new_id, bamfile]
            | map{ tuple(it[1], it[2]) } \
                //remove original id
            | groupTuple(by: 0) \
                // group by new_id
            | merge_dedup_bams
        sort(publishdir_ch, merge_dedup_bams.out)
        get_low_coverage(publishdir_ch, sort.out)

    emit:
        sort.out // bam_file channel

}
