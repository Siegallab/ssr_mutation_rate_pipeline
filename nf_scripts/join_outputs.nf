/*
 * Various functions for working with vcf, bed, and bam files
 */


process reorder_vcf_samples {

    // output a vcf with samples reordered according to sample_list

    input:
    tuple val(sample_list_str), path(vcf)

    output:
    path(ordered_vcf)

    script:
    if ("${vcf.extension}" == "gz") {
        ordered_vcf = "ordered_${vcf.name}"
    }
    else{
        ordered_vcf = "ordered_${vcf.name}.gz"
    }
    """
    module load ${params.BCFTOOLS}
    bcftools view -Oz -o $ordered_vcf --force-samples -s $sample_list_str $vcf
    """

}

process combine_vcf_regions {

    // combine freebayes outputs across regions
    // vcf files must already have samples in same order

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if "${publish_dir}" is empty string (interpreted as false)

    input:
    val publish_dir
    val file_prefix
    path(merge_input_file_list)

    output:
    path(out_vcf)

    script:
    out_vcf = "${file_prefix}_calls_comb.vcf.gz"
    concat_input = merge_input_file_list.join(' ')
    /*
     * need to include -a option to allow out-of-order concatenation 
     * of chromosome portions; -D probably not necessary
     * also not sure whether 'sort' necessary but probably can't hurt
     */
    """
    module load ${params.BCFTOOLS}
    echo $concat_input | xargs -n1 tabix -p vcf
    bcftools concat -Oz -o unsorted.vcf.gz -a -D $concat_input
    bcftools sort -Oz -o $out_vcf unsorted.vcf.gz
    """
}

process subset_vcf_by_reg {

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if "${publish_dir}" is empty string (interpreted as false)

    // subset vcf by region

    input:
    tuple val(region), path(vcf)

    output:
    tuple val(region), file(region_vcf)

    script:
    if ("$vcf" == params.fake_filename){
        publish_dir = ''
        region_vcf = params.fake_filename
        """
        if ! [ -f $region_vcf ]; then
            touch $region_vcf
        fi
        """
    }
    else{
        publish_dir = params.vcf_join_dir
        region_mod = region.replaceAll(/:/,'.') // change region call for filenames
        if ("${vcf.extension}" == "gz") {
            zip_command = ""
            zipped_vcf = vcf
            region_vcf = "filt_${region_mod}_${vcf}"
        }
        else {
            zip_command = "bgzip $vcf"
            zipped_vcf = "$vcf.gz"
            region_vcf = "filt_${region_mod}_${vcf}.gz"
        }
        """
        module load ${params.BCFTOOLS}
        $zip_command
        tabix -p vcf $zipped_vcf
        tabix -h $zipped_vcf $region | bgzip > $region_vcf
        """
    }
}

process subset_bam_by_reg {

    // subset bam by region

    input:
    tuple val(region), path(bam)

    output:
    tuple val(region), path(region_bam)

    script:
    region_mod = region.replaceAll(/:/,'.') // change region call for filenames
    region_bam = "chr_${region_mod}_${bam}"
    """
    module load ${params.SAMTOOLS}
    samtools index $bam
    samtools view -h -b $bam $region > $region_bam
    """

}

process create_target_bed {

    // remove params.excluded_regions from consideration

    input:
    val(region)

    output:
    tuple val(region), path(subset_bed)

    script:
    region_mod = region.replaceAll(/:/,'.') // change region call for filenames
    temp_bed_file = "region_${region_mod}_full.bed"
    subset_bed = "chr_${region_mod}.bed"
    """
    # create a bed file containing region
    # if region is full chromosome, start at 0, end at chromosome length
    if [[ $region == *\":\"* ]]; then
        chrom=\$(echo $region | sed 's/:.*//')
        start=\$(echo $region | sed 's/.*:\\([0-9]*\\).*/\\1/')
        if [[ $region == *\"..\"* ]]; then
            end=\$(echo $region | sed 's/.*:[0-9]*\\.\\.\\([0-9]*\\).*/\\1/')
        else
            end=\$(echo $region | sed 's/.*:[0-9]*-\\([0-9]*\\).*/\\1/')
        fi
    else
        chrom=$region
        start=0
        if [ ! -f ${params.ref}.fai ]; then
            module load ${params.SAMTOOLS}
            samtools faidx ${params.ref}
        fi
        end=\$(less ${params.ref}.fai | grep -P \"^\${chrom}\\t\" | awk '{print \$2}')
    fi
    echo -e \${chrom}\\\\t\${start}\\\\t\${end} > $temp_bed_file
    module load ${params.BEDTOOLS}
    bedtools subtract -a $temp_bed_file -b ${params.excluded_regions} > $subset_bed
    """

}

process normalize_vcf {

    /*
     * left-align, decompose, and normalize indels using bcftools
     * don't split multiple alt alleles into their own records
     * keep periods (not asterisks) for missing genotypes
     * avoids problems associated with vcflib
     */

    input:
    tuple path(in_vcf), val(split_alleles)

    output:
    path(out_vcf)

    script:
    out_vcf = "norm_${in_vcf.simpleName}.vcf.gz"
    if (split_alleles){
        m_option="-m-any"
    }
    else{
        m_option="-m+any"
    }
    """
    module load ${params.BCFTOOLS}
    # for normalization, need to replace any . in alt allele column with * 
    # per vcf format (doesn't seem to matter if we also replace records)
    bcftools view $in_vcf | \
        awk '{FS=OFS=\"\\t\"} {if (\$5 != \".\") gsub(/\\./,\"*\",\$5); print \$0}' | \
        bcftools norm --atomize --atom-overlaps . $m_option -f ${params.ref} -Oz -o ${out_vcf} -
    """
}

process filter_by_type {

    // filter vcf based on allele_type

    input:
    tuple path(in_vcf), val(allele_type)

    output:
    path(out_vcf)

    script:
    out_vcf = "filt_${in_vcf.simpleName}.vcf.gz"
    """
    module load ${params.BCFTOOLS}
    bcftools view --types $allele_type $in_vcf -Oz -o $out_vcf
    """
}

process vcf_to_csv {

    /*
     * converts vcf to per-sample csv (with samples that have missing
     * genotypes deleted), including 
     * CHROM,start,end,strain,GT,DP,allele,type,GL_diff
     */

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val(publish_dir)
    path(vcf) // list of input paths

    output:
    path(out_csv_zipped)

    script:
    out_csv = "${vcf.simpleName}.csv"
    out_csv_zipped = "${out_csv}.gz"
    """
    module load ${params.BCFTOOLS}
    sh ${params.shell_script_dir}/vcf_to_csv.sh $vcf $out_csv
    gzip $out_csv
    """
}

workflow combine_vcfs{

    // sort vcfs by sample and concatenate the sites

    take:
        publish_dir
        file_prefix
        sample_str_ch
            // channel with string of space-separated sample names
        vcf_ch

    main:
        reordered_vcf_ch = sample_str_ch \
            | combine(vcf_ch) \
            | reorder_vcf_samples \
            | collect
        combine_vcf_regions(publish_dir, file_prefix, reordered_vcf_ch)

    emit:
        combine_vcf_regions.out

}
