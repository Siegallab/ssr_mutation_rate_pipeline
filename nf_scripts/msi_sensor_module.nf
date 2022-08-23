/*
 * Perform analysis of SSR mutations using MSIsensor pro
 */

process create_default_msi_ssr_file {

    // creates ssr file that can be used for msi
    // uses default msisensor creator

    output:
    path(unfilt_ssr_file)

    script:
    unfilt_ssr_file = "msi_ssr_file.bed"
    """
    ${params.msi_sensor_pro} scan \
        -d ${params.ref} \
        -o $unfilt_ssr_file \
        -l ${params.min_homopolymer_rep} \
        -m 50 \
        -r ${params.min_nonhomopolymer_rep} \
        -s ${params.max_ssr_preprocess_period_size}
    """

}

process remove_excluded_msi {

    // remove excluded regions from msisensor input file

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val publish_dir
    path(unfilt_msi_bed)

    output:
    path(filt_msi_bed)

    script:
    prefilt_bed = "prefilt.bed"
    filt_msi_bed = "filt_${unfilt_msi_bed}"
    """
    # change heading (into new file)
    head -1 $unfilt_msi_bed | \
        awk '{FS=OFS=\"\\t\"} {print \$1,\$2,\"end\",\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' > \
        $prefilt_bed
    # exclude header, add end column and run bedtools to remove
    # params.excluded_regions from bed
    module load ${params.BEDTOOLS}
    tail -n +2 $unfilt_msi_bed | \
        awk '{FS=OFS=\"\\t\"} {end = \$2+\$3*\$5; print \$1,\$2,end,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' | \
        bedtools subtract -a - -b ${params.excluded_regions} >> $prefilt_bed
    # remove 'end' column
    cut -d\$'\\t' -f-2,4- $prefilt_bed > $filt_msi_bed
    """

}

process filter_strain_csv {

    // filter strain csv to exclude ancestral strains
    // remove all columns except strain and anc_strain, reverse their order
    // this can be done with channel manipilation but involves map,
    // which makes the process restart at every rerun

    input:
    path(strain_info_csv)

    output:
    path(strain_info_mod)

    script:
    strain_info_mod = "mod_$strain_info_csv"
    """
    awk 'BEGIN{FS=OFS=\",\"}{
        if (NR == 1) {
            for (i = 1; i <= NF; i++){
                if (\$i == \"strain\") {
                    strain_col = i
                }else if (\$i == \"anc_strain\") {
                    anc_strain_col = i
                }
            }
        }
        if (\$anc_strain_col != \$strain_col){
            print \$anc_strain_col,\$strain_col
        }
        }' $strain_info_csv > $strain_info_mod
    """

}

process run_msi_sensor {

    // run msi_sensor on refernce-sample pair

    input:
    tuple val(test_seq_id), val(ref_seq_id), path(ref_bam), path(test_bam), path(ssr_input)

    output:
    tuple val('summary'), path(out_summary)
    tuple val('somatic'), path(out_somatic)
    tuple val('germline'), path(out_germline)

    script:
    summary = "${test_seq_id}"
    somatic = "${test_seq_id}_somatic"
    germline = "${test_seq_id}_germline"
    out_summary = "${summary}.tsv"
    out_somatic = "${somatic}.tsv"
    out_germline = "${germline}.tsv"
    """
    module load ${params.SAMTOOLS}
    samtools index $ref_bam
    samtools index $test_bam
    ${params.msi_sensor_pro} msi \
        -d $ssr_input \
        -n $ref_bam \
        -t $test_bam \
        -o $summary \
        -c 15 \
        -l ${params.min_homopolymer_rep} \
        -p ${params.min_homopolymer_rep} \
        -m 120 \
        -q ${params.min_nonhomopolymer_rep} \
        -s ${params.min_nonhomopolymer_rep} \
        -w 120 \
        -f 0.1
    # add strain info
    awk '{FS=OFS=\"\\t\"}{if (NR==1) {print \$0,\"strain\"} else {print \$0,\"$test_seq_id\"}}' $summary > $out_summary
    awk '{FS=OFS=\"\\t\"}{if (NR==1) {print \$0,\"strain\"} else {print \$0,\"$test_seq_id\"}}' $somatic > $out_somatic
    awk '{FS=OFS=\"\\t\"}{if (NR==1) {print \$0,\"strain\"} else {print \$0,\"$test_seq_id\"}}' $germline > $out_germline
    """

}

process concat_files {

    // concatenates output files of a specific category

    // creates ssr file that can be used for msi

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val publish_dir
    tuple val(category), path(files)

    output:
    path(cat_file)

    script:
    first_file = files[0]
    file_list = files.join(' ')
    cat_file = "msi_out_${category}_comb.tsv"
    """
    # write heading line from any file
    head -1 $first_file > $cat_file
    # write non-heading lines from all files
    tail -n +2 -q $file_list >> $cat_file
    """

}

workflow run_msi {

    // perform post-alignment duplicate removal and merging based on id

    take:
        publishdir_ch
            // Value channel containing directory for alignments;
            // if empty, alignments not saved
        sort_out_ch // [id, bam]
        ssr_file_path

    main:
        msi_ssr_unfilt_ch = create_default_msi_ssr_file()
        msi_ssr_path = remove_excluded_msi(publishdir_ch, msi_ssr_unfilt_ch)

        Channel.fromPath("${params.strain_seq_info_dir}/${params.strain_info_file}") \
            | filter_strain_csv \
            | splitCsv(header: false, skip: 1) \
                // [ref_id, ma_id]
            | set { ref_dict_ch }
        msi_out_ch = \
            ref_dict_ch \
            | combine(sort_out_ch, by:0) \
                // [ref_id, ma_id, ref_bam]
            | map{ [it[1], it[0], it[2]]} \
            | combine(sort_out_ch, by:0) \
                // [ma_id, ref_id, ref_bam, ma_bam]
            | combine(msi_ssr_path) \
                // [ma_id, ref_id, ref_bam, ma_bam, msi_ssr_path]
            | run_msi_sensor \
            | mix \
            | groupTuple

        concat_files(publishdir_ch, msi_out_ch)

    emit:
        concat_files.out // list of germline, somatic, and summary output files

}