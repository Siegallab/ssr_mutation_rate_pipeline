#!/bin/bash

in_vcf=$1
out_csv=$2
normalize_indel_alleles=${3:-F}


# need to check that vcf isn't empty
if [ -s $in_vcf ]; then
    # normalizing indels (bcftools 1.14+) allows bcftools query to work with
    # long lists of long alleles and extract %TGT value from them without 
    # running into error, but is probably not necessary in most cases
    if [[ $normalize_indel_alleles == T ]]; then
        # for normalization, need to replace any . in alt allele column with * 
        # per vcf format (doesn't seem to matter if we also replace records)
        bcftools view $in_vcf | \
            awk '{FS=OFS="\t"} {gsub(/\./,"*",$5); print $0}' | \
            bcftools norm -a - -f $ref -o norm_${in_vcf}

        # if doing normalization, need to use normalized vcf as input
        data_extraction_input=norm_${in_vcf}
    else
        data_extraction_input=${in_vcf}
    fi

    # extract GL values (inherently comma-separated) for every sample
    bcftools query -f '[%GL\n]' ${data_extraction_input} > GL.csv
    # set GL_diff to max of GL values, excluding a single zero 
    # (corresponding to ML allele)
    # This ensures that if top two alleles have same likelihood, GL_diff 
    # is reported as 0
    awk 'BEGIN{FS=OFS=","}
        {
            curr_max = -1000000
            for (j=1; j <= NF; j++)
                if( $j > curr_max )
                    curr_max = $j
            GL_diff = -1000000
            max_traversed = 0
            for (i = 1; i <= NF; ++i){
                curr_diff = $i - curr_max
                if ($i == curr_max && max_traversed == 0)
                    {max_traversed = 1; max_i = i}
                else if (curr_diff > GL_diff)
                    GL_diff = curr_diff
                }
            print GL_diff
        }' GL.csv > GL_diff.csv

    # extract AD values (inherently comma-separated) corresponding to 
    # GT for every sample
    # only fill in AD column if GT column doesn't contain | (i.e. haploids)
    bcftools query -f '[%GT\t%AD\n]' ${data_extraction_input} | \
        awk 'BEGIN{FS="\t"; OFS=","} {split($2,a,/,/); if ( !($1 ~ /"|"/) && !($1 ~ /"\/"/) ) { print a[$1+1] } }' >> AD.csv

    # only save if GT is non-missing
    # some vcfs don't have 'TYPE', just don't include it if that's the case
    bcftools query -f '[%CHROM\t%POS\t%END\t%SAMPLE\t%GT\t%DP\t%TGT\t%INFO/TYPE\t%REF\t%GQ\n]' ${data_extraction_input} > temp_vcf_data_pre_awk.csv
    if [ $? -eq 0 ]; then
        awk 'BEGIN{FS="\t"; OFS=","} {split($8,a,","); gsub(/,/,";",$8); if ($5==0) {print $1,$2,$3,$4,$5,$6,$7,"",$8,$9,$10} else {print $1,$2,$3,$4,$5,$6,$7,a[$5],$8,$9,$10}}' temp_vcf_data_pre_awk.csv > temp_vcf_data.csv
        echo "CHROM,start,end,strain,GT,DP,allele,allele_type,type,ref_allele,GQ,GL_diff,AD" > $out_csv
    else
        bcftools query -f '[%CHROM\t%POS\t%END\t%SAMPLE\t%GT\t%DP\t%TGT\t%REF%GQ\n]' ${data_extraction_input} > temp_vcf_data_pre_awk.csv
        awk 'BEGIN{FS="\t"; OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' temp_vcf_data_pre_awk.csv > temp_vcf_data.csv
        echo "CHROM,start,end,strain,GT,DP,allele,ref_allele,GQ,GL_diff,AD" > $out_csv
    fi
    paste -d',' temp_vcf_data.csv GL_diff.csv AD.csv | \
        awk -F"," '$5 != "\." { print $0 }' >> $out_csv
    rm temp_vcf_data.csv
    rm temp_vcf_data_pre_awk.csv
    rm GL_diff.csv
    rm GL.csv
    rm AD.csv
else
    touch $out_csv
fi

