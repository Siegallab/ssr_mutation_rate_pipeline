/*
 * Perform post-analysis of called mutations in R
 */

process postanalyze_call_results {

    // output .Rdata file with results of R analysis

    publishDir path: "${params.out_dir}/${params.final_calls_dir}", mode:'copy', overwrite: true

    input:
    val(r_out_prefix)
    path(ssr_bed_file)
    path(call_csv)

    output:
    path(r_out)

    script:
    r_out = "${r_out_prefix}_analyzed.Rdata"
    """
    module load ${params.R}
    Rscript --vanilla ${params.r_script_dir}/seq_postanalyze.R \
        ${params.r_script_dir}/supp_calling_functions.R \
        $call_csv \
        ${params.ref} \
        ${params.excluded_regions} \
        ${params.strain_seq_info_dir}/${params.strain_info_file} \
        $ssr_bed_file \
        ${params.ploidy} \
        ${params.max_BQ} \
        ${params.rel_call_r_fun} \
        ${params.slide_window_size_str} \
        ${params.ssr_grouping_cat} \
        ${params.additional_filter_exp} \
        ${params.tel_cen_ltr_exclusion_dist} \
        $r_out
    """

}