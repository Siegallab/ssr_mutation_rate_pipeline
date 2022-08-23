/*
 * Create and clean bedfile containing locations of SSR loci
 *
 * Python portions require installation of biopython, pyvcf, numpy
 */


include { combine_vcfs } from './join_outputs.nf'

process split_fasta {

    // split fasta file by chromosome
    input:
    path(ref)

    output:
    path("MA-chrom*.fasta")

    script:
    """ 
    awk -F \"|\" '/^>/ \
        {close(F); ID=\$1; gsub(\"^>\", \"\", ID); gsub(\" dna.*\", \"\", ID); F=\"MA-chrom\"ID\".fasta\"} \
        {gsub(\" dna.*\$\", \"\"); print >> F}' \"${ref}\"
    """
}

process create_ssr_pattern {

    // create pattern file for ssr search

    output:
    path(ssr_pattern_file)

    script:
    ssr_pattern_file = "ssr_pattern.fa"
    """
    module load ${params.python}
    python ${params.py_script_dir}/write_ssr_pattern_fasta.py ${ssr_pattern_file}
    """

}

process run_trf {

    /*
     * identify perfect ssr loci <= 10 bp & with motif length <= 3
     * by direct search on single-chrom fasta file
     */
    /* combine with trf output */
    /* perform filtration on results and split overlapping loci */

    input:
    path(fasta)

    output:
    path(out_file)

    script:
    // Default parameters to use with Tandem Repeat Finder
    // These are the recommended values
    MATCH_WT=params.match_score
    MISMATCH_PEN=params.mismatch_penalty
    INDEL_PEN=params.indel_penalty
    P_MATCH=80
    P_INDEL=10
    MAX_PERIOD=500
    chrom = (fasta =~ /MA-chrom(.*)\.fasta/)[0][1]
    out_file = "trf_out_fixed_${chrom}.txt"
    """
    module load ${params.TRF}
    trf409.linux64 \
        $fasta $MATCH_WT $MISMATCH_PEN $INDEL_PEN $P_MATCH $P_INDEL \
        ${params.min_trf_alignment_score} $MAX_PERIOD \
        -h -d -l 6 -ngs > trf_out_init.txt
    # filter out repeats with a period longer than max_ssr_preprocess_period_size and fix a 
    # few issues with incorrect TRF entries
    module load ${params.python}
    python ${params.hipstr_py_script_dir}/fix_trf_output_mod.py trf_out_init.txt ${out_file} ${params.max_ssr_preprocess_period_size}
    """
}

process id_short_ssr_loci {

    /*
     * identify perfect ssr loci <= 10 bp & with motif length <= 3
     * by direct search on single-chrom fasta file
     */

    /*
     * nb: here, we don't filter based on max period size, since we 
     * assume it's longer than 3 (the max period size created by 
     * this process)
     */

    input:
    tuple path(fasta), path(ssr_pattern_file)

    output:
    path(out_file)

    script:
    chrom = (fasta =~ /MA-chrom(.*)\.fasta/)[0][1]
    out_file = "short_ssr_loci_${chrom}.txt"
    """
    # start file with chrom, as in TRF
    echo '@${chrom}' > ${out_file}
    # need to convert start and end to same indexing as trf
    # unneeded columns will be filled in with dashes
    ${params.seqkit} locate -Pri -f $ssr_pattern_file $fasta --bed | \
        awk 'BEGIN{FS=\"\\t\"; OFS=\" \"}
            {   start=\$2+2;
                end=\$3-1;
                split(\$4,info,/-/);
                motif = info[1];
                motif_len = info[2];
                num_copies = info[3];
                score = info[4];
                ssr_seq = info[5];
                {
                    print start,end,motif_len,num_copies,motif_len,100,0,score,\"-\",\"-\",\"-\",\"-\",\"-\",motif,ssr_seq,\"-\",\"-\"
                }
            }' >> $out_file
    """
}

process filter_repeat_bed {

    /*
     * combine repeat bed files from multiple chromosomes, reformat the 
     * TRF entries and filter to only include repeats with a 
     * sufficiently high score and shorter than params.max_ssr_length
     */

    input:
    tuple path(file_list), val(repeat_filetype) // list of input paths

    output:
    path(joined_bed)

    script:
    joined_bed = "${repeat_filetype}_joined_filtered.bed"
    file_list_str = file_list.join(',')
    """
    # reformat the TRF entries and filter to only include repeats 
    # with a sufficiently high score:
    module load ${params.python}
    python ${params.hipstr_py_script_dir}/trf_parser_mod.py \
        $file_list_str \
        ${params.max_ssr_preprocess_period_size} \
        ${params.min_homopolymer_rep} \
        ${params.min_nonhomopolymer_rep} > $joined_bed
    """


}

process combine_bed {

    /*
     * run additional hipst-suggested filtration, combine results into 
     * single bed file
     */

    publishDir path: "${params.out_dir}/${publish_dir}", mode:'copy', overwrite: true, enabled: "${publish_dir}"
        // doesn't publish if ${publish_dir} is empty string (interpreted as false)

    input:
    val(publish_dir)
    val(merge_neighbors)
    tuple path(trf_bed), path(short_ssr_bed)

    output:
    path(out_bed)

    script:
    out_bed="ssr_preprocess_reference_${params.genome_name}_sitemerge_${merge_neighbors}.bed"
    filtered_repeat_file = "filtered_repeats.bed"
    overlap_pass = "overlap_pass.bed"
    overlap_fail = "overlap_fail.bed"
    final_merged_file = "merge_pass.bed"
    if (merge_neighbors) {
        merging_string = ""
        merge_site_code = \
            """
            module load ${params.python}
            python ${params.hipstr_py_script_dir}/analyze_overlaps.py filtered_repeats.sorted.bed $overlap_pass $overlap_fail
            ########################################################################
            # cleanup overlaps just in case, to make sure end is not smaller than beginning
            awk 'BEGIN{F=OFS=\"\\t\"}  {if (\$3 >= \$2) {print \$0}}' $overlap_pass > temp_$overlap_pass
            awk 'BEGIN{F=OFS=\"\\t\"}  {if (\$3 >= \$2) {print \$0}}' $overlap_fail > temp_$overlap_fail
            mv temp_$overlap_pass $overlap_pass
            mv temp_$overlap_fail $overlap_fail
            # remove any entries within 10bp of a failed merge region:
            bedtools window -w 10 -a $overlap_pass -b $overlap_fail -v \
                > pass.r2
            # To minimize the effects of nearby STRs on genotyping errors,
            # extract entries that aren't within 10bp of another entry or are
            # within 10bp of one or more entries that all share the same period
            bedtools merge -i pass.r2 -c 4,6 -o collapse -d 10 \
                | grep -v \",\" > $final_merged_file
            bedtools merge \
                -i pass.r2 \
                -c 4,4,4,6 \
                -o collapse,count_distinct,distinct,collapse \
                -d 10 \
                | grep \",\" \
                | awk '\$5 == 1' \
                | awk -v OFS=\"\\t\" '{print \$1, \$2, \$3, \$6, \$7}' \
                | sed \"s/,/\\//g\" \
                >> $final_merged_file
            """
    }
    else {
        merge_site_code = \
            """
            module load ${params.python}
            python ${params.hipstr_py_script_dir}/analyze_overlaps_split_motifs.py \
            filtered_repeats.sorted.bed \
            $overlap_pass $overlap_fail \
            ${params.min_homopolymer_rep} \
            ${params.min_nonhomopolymer_rep} \
            ${params.match_score} \
            ${params.mismatch_penalty} \
            ${params.indel_penalty} \
            ${params.min_score_per_bp}
           cat $overlap_fail $overlap_pass | \
                awk 'BEGIN{F=OFS=\"\\t\"}  {if (\$3 >= \$2) { print \$1,\$2,\$3,\$4,\$6 } }' | \
                bedtools sort -i stdin > $final_merged_file
            """
        merging_string = "un"
    }
    """
    cat $short_ssr_bed $trf_bed > $filtered_repeat_file
    ########################################################################
    module load ${params.BEDTOOLS}
    bedtools sort -i $filtered_repeat_file | \
        uniq > \
        filtered_repeats.sorted.bed
    ########################################################################
    # merge overlapping STRs into single entries or split STRs to not overlap
    $merge_site_code
    # Lastly, construct the final reference
    # the script below assumes there's no tabs? (I think) in the bed 
    # file labels
    cat $final_merged_file \
        | bedtools sort \
        | awk -v merging_str=${merging_string} \
        'BEGIN{OFS=\"\\t\"} {print \$1, \$2, \$3, \$4, (\$3-\$2+1)/\$4, \"YEAST_STR_MOD_\"merging_str\"merged_\"NR, \$5}' > \
        $out_bed
    """

}

process remove_called_repeats {

    // remove any loci that overlap with loci in called_vcf from bedfile bed_in

    input:
    val(bed_zero_indexed)
    path(called_vcf)
    path(bed_in)

    output:
    path(bed_out)

    script:
    bed_out="call_subtracted_repeats.bed"
    bed_from_vcf="bed_from_vcf.bed"
    if (bed_zero_indexed){
        vcf_idx_correction='1'
    }
    else{
        vcf_idx_correction='0'
    }
    """
    if [ ! -f ${params.ref}.fai ]; then
        module load ${params.SAMTOOLS}
        samtools faidx ${params.ref}
        module purge
    fi
    # bedtools has different convention than vcf
    # bed files are 0-indexed and the end position is actually one after the end
    # vcf files are 1-indexed and the end is part of the sequence
    # therefore, first convert .vcf to .bed, adding 1 to end and
    # converting to indexing corresponding to that of $bed_in
    # Merge is tricky: we want to do two things:
    #   1.  Remove any SSRs where the ref allele intersects with the 
    #       SSR position
    #   2.  Remove any SSRs where the ref allele is the nucleotide 
    #       before the SSR start position, and the alt allele type
    #       includes an insertion
    # NB: this may include non-SSR-motif insertions, but looking for 
    # only SSR motif-related insertions is probably too complicated
    module load ${params.BEDTOOLS}
    module load ${params.BCFTOOLS}
    bcftools query -f '[%CHROM\\t%POS\\t%END\\t%REF\\t%TGT\\n]' $called_vcf | \
        awk 'BEGIN{FS=OFS=\"\\t\"} {print \$1,\$2-$vcf_idx_correction,\$3-$vcf_idx_correction+1,\$4,\$5}' > \
        $bed_from_vcf
    # get SSR loci that don't intersect with called loci in vcf
    bedtools subtract -A -a $bed_in -b $bed_from_vcf > internal_nonmatch_$bed_out
    # Among those SSR loci that don't intersect with called loci in vcf,
    #   1.  convert them to only include nucleotide before SSR;
    #       keep motif as column
    #   2.  intersect with $bed_from_vcf
    #   3.  keep only those rows in which one of the called alleles is 
    #       an insertion
    awk 'BEGIN{FS=OFS=\"\\t\"} {
            if (\$2 == 0) {start = 0; end = 1} else {start = \$2-1; end = \$2}
            {print \$1,start,end,\$2,\$3,\$4,\$5,\$6,\$7}
        }' internal_nonmatch_$bed_out | \
        bedtools intersect -wa -wb -a stdin -b $bed_from_vcf | \
        awk 'BEGIN{FS=OFS=\"\\t\"} {if ( length(\$14) > length(\$13) ) {print \$1,\$4,\$5}}' | \
        bedtools sort -i stdin | uniq > start_overlap_$bed_out
    bedtools subtract -A -a internal_nonmatch_$bed_out -b start_overlap_$bed_out | \
        bedtools sort -i stdin -faidx ${params.ref}.fai \
        > $bed_out
    """

}

process bed_to_vcf {

    /*
     * convert (trf-generated, postprocessed) bed file to vcf file with 
     * reference alles (and missing alt alleles)
     */

    input:
    val(zero_indexed)
    tuple path(ref_fasta), path(bed_file)

    output:
    path(out_vcf)

    script:
    if (zero_indexed){
        bed_idx_correction='0'
    }
    else{
        bed_idx_correction='1'
    }
    out_vcf = "ssr_preprocess_reference_${params.genome_name}_bed.vcf.gz"
    intermed_file = "temp.bed"
    """
    module load ${params.BEDTOOLS}
    module load ${params.BCFTOOLS}
    # get only sequence of nucleotide before repeat followed by first 
    # motif copy's worth of nucleotides; this will make number of 
    # nucleotides used in call more similar to what freebayes uses when 
    # making a real call
    # need to make adjustment for case in which motif starts at first 
    # nucleotide in chromosome
    awk  -F'\\t' \
        'BEGIN{OFS = FS} {
            if (\$2 - $bed_idx_correction >= 1)
            {mod_start = \$2 - $bed_idx_correction - 1}
            else
            {mod_start = \$2 - $bed_idx_correction}
            {print \$1, mod_start, mod_start+1+\$4}
        }' ${bed_file} | \
        bedtools getfasta -fi ${ref_fasta} -bed - -tab | \
        awk -F'\\t' '{print \$2}' > ref_seqs.txt
    # combine reference sequence with the rest of bed file into vcf
    # alt alleles etc set to '.'
    paste -d'\\t' ref_seqs.txt ${bed_file} | \
        awk -F'\\t' 'BEGIN{OFS = FS} {
            ori_bed_start = \$3 - $bed_idx_correction;
            ori_vcf_start = ori_bed_start+1;
            if (ori_vcf_start > 1){mod_vcf_start = ori_vcf_start-1}
            else {mod_vcf_start = ori_vcf_start}
            {print \$2,mod_vcf_start,\$7,\$1,\".\",\".\",\".\",\".\"}
            }' > \
        ${intermed_file}
    echo -e \"#CHR\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\" | \
        cat - ${intermed_file} | bgzip > ${out_vcf}
    """

}

workflow create_bed {

    take:
        bedfile_publishdir_ch
        merge_neighbors_ch

    main:
        /*
         * merge_neighbors_ch determines whether bed files should be 
         * created according to hipstr specs (merging nearby repeats) 
         * or not (restricting trf output to only top repeat at a 
         * locus, evenly splitting all repeats with any overlap, no 
         * merging of neighbors)
         */

        ref_ch = channel.fromPath(params.ref)
        ssr_pattern_ch = create_ssr_pattern()
        split_fasta_ch = \
            ref_ch \
            | split_fasta \
            | flatten
        trf_name_ch = channel.value('trf')
        ssr_name_ch = channel.value('short_ssr')
        trf_out_ch = \
            split_fasta_ch \
            | run_trf \
            | combine(trf_name_ch) \
            | groupTuple(by:1)
        short_ssr_out_ch = \
            split_fasta_ch \
            | combine(ssr_pattern_ch) \
            | id_short_ssr_loci \
            | combine(ssr_name_ch) \
            | groupTuple(by:1)
        filtered_repeat_ch = \
            trf_out_ch \
            | concat(short_ssr_out_ch) \
            | filter_repeat_bed \
            | collect
        combine_bed(bedfile_publishdir_ch, merge_neighbors_ch, filtered_repeat_ch)

    emit:
        combine_bed.out

}

workflow ssr_preprocess {

    take:
        bedfile_publishdir_ch
        called_vcf_ch

    main:
        merge_neighbors_ch = channel.value(false)
        bed_file_ch = create_bed(bedfile_publishdir_ch, merge_neighbors_ch)
        zero_idx_ch = channel.value(false)
        uncalled_bed_file_ch = \
            remove_called_repeats(zero_idx_ch, called_vcf_ch, bed_file_ch)
        ref_ch = channel.fromPath(params.ref)
        vcf_to_convert_ch = \
            ref_ch \
            | combine(uncalled_bed_file_ch)
        bed_to_vcf(zero_idx_ch, vcf_to_convert_ch)
        out_ch = bed_file_ch.concat(bed_to_vcf.out)

    emit:
        out_ch

}