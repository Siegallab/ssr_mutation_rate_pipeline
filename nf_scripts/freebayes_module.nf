/*
 * call mutations using freebayes
 */

include { subset_vcf_by_reg; subset_bam_by_reg; create_target_bed; combine_vcfs } from './join_outputs.nf'

process run_freebayes {

    /* 
     * run freebayes across samples, in $region, with input files 
     * specified by $bamlist
     */
    // BAMs should be sorted

    input:
    val filter
    tuple val(region), path(bamlist), path(vcf_variant_input), path(target_bed)

    output:
    path(out_vcf)

    script:
    region_mod = region.replaceAll(/:/,'.') // change region call for filenames
    out_vcf = "freebayes_calls_chrom_${region_mod}.vcf.gz"
    /*
     * only include --variant_input option if vcf_variant_input isn't
     * empty string
     */
    if ("${vcf_variant_input.name}" == params.fake_filename){
        variant_input_line = ''
        vcf_subset_line = ''
        tabix_line = ''
    }
    else {
        vcf_file_extension = vcf_variant_input.getExtension()
        if ("${vcf_variant_input.extension}" == "gz") {
            zipped_vcf = vcf_variant_input
            tabix_line = "tabix -p vcf $zipped_vcf"
        }
        else {
            zipped_vcf = "${vcf_variant_input}.gz"
            tabix_line = "bgzip $vcf_variant_input; tabix -p vcf $zipped_vcf"
        }
        variant_input_line = "--variant-input $zipped_vcf --only-use-input-alleles --min-alternate-count 0 --min-alternate-fraction 0 --min-coverage 0 "
    }
    bam_input = bamlist.join(' --bam ')
    bam_list_str_xargs = bamlist.join(' ')
    if (filter == ""){
        post_fb_code = \
        """
        module purge
        module load ${params.SAMTOOLS}
        cat unfiltered.vcf | bgzip > $out_vcf
        """
    }
    else{
        post_fb_code = \
        """
        module purge
        module load ${params.VCFLIB}
        module load ${params.SAMTOOLS}
        vcffilter -f \"${filter}\" unfiltered.vcf | bgzip > $out_vcf
        """
    }
    """
    # recreate index files locally
    module load ${params.SAMTOOLS}
    $tabix_line
    #echo -n $bam_list_str_xargs | xargs -d \" \" -I@ samtools sort -O BAM -o @ @
    echo $bam_list_str_xargs | xargs -n1 samtools index
    if [ ! -f ${params.ref}.fai ]
    then
        samtools faidx ${params.ref}
    fi
    # run freebayes
    ${params.FREEBAYES_BINARY} \
        --bam ${bam_input} \
        --fasta-reference ${params.ref} \
        --target ${target_bed} \
        --ploidy ${params.ploidy} \
        --report-genotype-likelihood-max \
        --min-mapping-quality 1 \
        --genotype-qualities $variant_input_line \
        --vcf unfiltered.vcf
    # apply filter if necessary and zip
    $post_fb_code
    """
}

workflow freebayes {

    take:
        publishdir_ch
            // Value channel containing directory for freebayes output
            // if empty, no reference copied
        prefix_ch
            // value channel with prefix to add before fb vcf output file name
        bam_ch // bamfile channel (should already by sorted)
        vcf_variant_input_ch // value channel with path to vcf file
        filter_ch // value channel with filter to apply to freebayes results
        sample_list_ch // value channel with string of comma-separated samples
    main:
        /*
         * We will run freebayes on individual chromosomes to speed up
         * processing time, and then merge the results
         *
         * Note we remove 'Mito' from the list of choromosomes because 
         * FreeBayes fails at analyzing it
         * Also, analysis of chroms XII takes forever, so split it into three
         * (excluding most of the rDNA region)
         */
        chrom_ch = channel.fromList( \
            params.region_list
            )
        // bams should already be sorted
        /*
         * There is seemingly a bug in FreeBayes that means that if a 
         * vcf variant input is supplied, bam files and vcf files need 
         * to be broken up by chromosome
         */
        vcf_variant_input_reg_ch = \
            chrom_ch \
            | combine(vcf_variant_input_ch) \
            | subset_vcf_by_reg
        bam_reg_ch = \
            chrom_ch \
            | combine(bam_ch) \
            | subset_bam_by_reg \
            | groupTuple(by: 0)
                // [region, [bam_list]]
        /* 
         * also create 'target' bed file that includes current region 
         * and excludes params.excluded_region_bed
         */
        target_reg_ch = \
            chrom_ch \
            | create_target_bed
        fb_in_ch = \
            bam_reg_ch \
            | combine(vcf_variant_input_reg_ch, by: 0) \
            | combine(target_reg_ch, by: 0)
                // [region, [bam_list], vcf_variant_input, target_bed]
        run_freebayes(filter_ch, fb_in_ch) \
                // vcf_file
            | set { fb_out_ch }
        combine_vcfs(publishdir_ch, prefix_ch, sample_list_ch, fb_out_ch)

    emit:
        combine_vcfs.out // vcf_file

}