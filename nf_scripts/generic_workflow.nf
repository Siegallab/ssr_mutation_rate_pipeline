/*  
 *
 * Pipeline to perform alignment, run freebayes, update reference 
 * genome using identified variants, re-run alignment, and run 
 * freebayes using reads aligned to updated reference genome
 *
 */

include { trim; idx_for_align; align_mem; postalign_process; subset_fastq } from './prep_align_module.nf'

include {run_msi} from './msi_sensor_module.nf'

include { freebayes; freebayes as re_bayes } from './freebayes_module.nf'

include { combine_vcfs; normalize_vcf; normalize_vcf as normalize_uncalled_vcf; filter_by_type; vcf_to_csv; vcf_to_csv as vcf_to_csv_ii } from './join_outputs.nf'

include { ssr_preprocess } from './ssr_preprocess_module.nf'

include { postanalyze_call_results } from './postanalysis.nf'

workflow align_call_join_recall {

    // Combine alignment, variant calling, joining variant calls, and re-analyzing variant calls

    take:
        input_file_ch

    main:
        // declare channels to be used for output directories
        align_dir_ch = channel.value(params.aligned_out)
        final_calls_dir_ch = channel.value(params.final_calls_dir)
        blank_ch = channel.value("")
        sample_list_ch = \
            channel.value(GroovyCollections.transpose(params.id_dict)[1].unique(false).join(','))
        // subset reads if necessary
        if (params.containsKey('subset_read_num'))
            trim_input_ch = input_file_ch | subset_fastq
        else
            trim_input_ch = input_file_ch
        // trim reads
        trim(trim_input_ch) | set {trim_out_ch}
        // Run align
        channel.value(params.ref) \
            | idx_for_align \
            | combine(trim_out_ch) \
            | align_mem \
            | set { init_align_out_ch }
        // postprocess alignments        
        sorted_out_ch = postalign_process(align_dir_ch, init_align_out_ch)
        postaligned_bam_ch = 
            sorted_out_ch \
            | toList \
            | transpose \
            | last \
            | flatten
        // make initial fb calls
        init_call_dir_ch = channel.value(params.initial_calls_dir)
        init_vcf_variant_input_ch = channel.fromPath(params.fake_filename)
        init_fb_filter_ch = channel.value("QUAL > 1")
        split_allele_init_ch = channel.value(false)
        init_fb_calls_out_ch = \
            freebayes(init_call_dir_ch, \
                channel.value("init_call_freebayes"), \
                postaligned_bam_ch, \
                init_vcf_variant_input_ch, \
                init_fb_filter_ch,
                sample_list_ch) \
            | combine(split_allele_init_ch) \
            | normalize_vcf
        // convert fb output to csv
        vcf_to_csv(init_call_dir_ch, init_fb_calls_out_ch)
        // Create vcf with reference alleles for uncalled, unmerged SSRs
        ssr_preprocess(final_calls_dir_ch, init_fb_calls_out_ch)
        fb_ssr_bed_ch = ssr_preprocess.out.first()
                genotyping_template_uncalled_vcf_ch = ssr_preprocess.out.last()
                // make msi calls
                run_msi(
                    channel.value(params.msisensor_out_dir),
                    sorted_out_ch,
                    fb_ssr_bed_ch
                    )
        /*
         * feed alignment + uncalled ssr template vcf into freebayes
         * 
         * freebayes takes ref-only calls in vcf here because it's acting as
         * a genotyper for all ssr loci that don't have a mutation detected
         * 
         * freebayes treats missing alt alleles in --only-use-input-alleles
         * by converting them to two possible alleles: a deletion of the 
         * first nucleotide, and a deletion of everything but the first 
         * nucleotide (in this code, this would be one motif)
         * 
         * it seems that the GL values for such alt alleles for loci 
         * where no alleles were called by the initial caller are 
         * identical; unclear whether this is real (i.e. both alleles 
         * are equally unlikely) or due to freebayes bug that causes it 
         * to treat alt alleles identically if one of the supplied alt 
         * alleles is unreasonable (no support)
         */
        split_allele_uncalled_ch = channel.value(true)
        uncalled_fb_calls_out_ch = \
            re_bayes(blank_ch, \
                channel.value("uncalled_freebayes"), \
                postaligned_bam_ch, \
                genotyping_template_uncalled_vcf_ch, \
                blank_ch, \
                sample_list_ch) \
            | combine(split_allele_uncalled_ch) \
            | normalize_uncalled_vcf \
            | combine(channel.value('indels')) \
            | filter_by_type
        /*
         * combine and sort outputs of two freebayes runs
         * here, we will use combine_vcf_regions to do simple 
         * bcftools concatenation
         */
        fb_calls_ch = \
            init_fb_calls_out_ch \
            | mix(uncalled_fb_calls_out_ch)
        combined_fb_calls_ch = \
            combine_vcfs(final_calls_dir_ch, channel.value("fb"), sample_list_ch, fb_calls_ch)
        // convert fb output to csv
        vcf_to_csv_ii(final_calls_dir_ch, combined_fb_calls_ch)
        // Run R-based postanalysis
        postanalyze_call_results('fb',fb_ssr_bed_ch,vcf_to_csv_ii.out)

}
