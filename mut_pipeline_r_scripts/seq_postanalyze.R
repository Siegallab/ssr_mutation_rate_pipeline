################################################################################
################################################################################
# script to analyze output of sequencing data
################################################################################
################################################################################
# accept external arguments
args = commandArgs(trailingOnly=TRUE)

script_file <- args[1]
data_csv <- args[2]
genome_file <- args[3]
excluded_regions_bedfile <- args[4]
strain_info_file <- args[5]
bed_file_full <- args[6]
ploidy <- as.integer(args[7])
QUAL <- as.integer(args[8])
rel_call_fun_name <- args[9]
slide_window_size_str <- args[10]
ssr_grouping_cat <- args[11]
additional_filter_exp <- args[12]
tel_cen_ltr_exclusion_dist <- as.numeric(args[13])
out_file <- args[14]

################################################################################
################################################################################

library(conflicted)
library(Matrix)
library(stringr)
library(seqinr)
library(tidyverse)

################################################################################
################################################################################
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rowRanges", "DelayedArray")
options(dplyr.summarise.inform = FALSE)
################################################################################
################################################################################
source(script_file)
rel_call_fun <- get(rel_call_fun_name)
################################################################################
# set some parameters
################################################################################
# number of bins to divide AT proportion approximation into
# (rounds to nearest (1-bins)th)
AT_prop_approx_bins = 3

# distance for mutations to be considered as making up a single 'superlocus'
non_SSR_superlocus_dist <- 50

# min number of ssr loci for a 'category' GL_diff_threshold to not be thrown out
# at below 26, poorly genotyped loci with low GL_diff thresholds get through
# resulting in clearly incorrect calls
# (however, it's possible those can be filtered out by removing loci mutated
# across multiple strains?)
min_ssr_locus_per_category <- 26

# convert QUAL from phred score to probability
p_qual <- 10^(-as.numeric(QUAL)/10)

# properties of SSR quantile testing
quant_list <- round(seq(0,1,by = 0.05),2)
slide_window_sizes <- as.numeric(unlist(strsplit(slide_window_size_str,',')))

# min length (in bp) of an SSR for it to be considered 'long'
long_ssr_min_length <- 8

################################################################################
# read in excluded regions, ssr regions, and strain info
################################################################################

# add 1 to excluded region start to offset bedfile 0-indexing
excluded_regions_df <-
  read_tsv(
    excluded_regions_bedfile,
    col_names = c('CHROM', 'start', 'end')
    ) %>%
  mutate(start = start + 1) %>%
  rowid_to_column("locus") # set arbitrary locus column

strain_df <-
  read_csv(strain_info_file)

strain_anc_df <-
  strain_df %>%
  filter(cell_divisions > 0)

# Load in .bed file
ssr_bed_df_full <- read_tsv(bed_file_full, col_names = F)
colnames(ssr_bed_df_full) <- c('CHROM','start','end','MOTIF_LEN','NUM_COPIES','locus','MOTIF')

################################################################################
# Load in data
################################################################################
curr_data <-
  read_csv(data_csv,
           na = c("",".","NA"),
           col_types = cols(allele_type = col_character())) %>%
  mutate(locus = paste0(CHROM,":",start,"-",end)) %>%
  group_by(locus) %>%
  mutate(
    unique_types = 
      ifelse(all(is.na(allele_type)),
             NA,
             paste(sort(unique(allele_type[!is.na(allele_type)])),collapse = ';')
      )
  )%>%
  ungroup()

pos_df <-
  curr_data %>%
  distinct(locus,CHROM,start,end,ref_allele) %>%
  mutate(long_allele = nchar(ref_allele) > 1)

################################################################################
# Identify SSRs, excluded regions
################################################################################
ssr_match_df_by_start <-
  recode_loci_overlap(
    pos_df,
    mutate(
      ssr_bed_df_full,
      start = pmax(start-1,0),
      end = start+1+MOTIF_LEN
    ),
    alt_chrom_key = list('Mito' = 17*10^7),
    bed_dist_expansion = 0
  ) %>%
  dplyr::rename(ssr = bed_locus, locus = pos_locus)

ssr_match_df_unneat <-
  recode_loci_overlap(
    pos_df,
    ssr_bed_df_full,
    alt_chrom_key = list('Mito' = 17*10^7),
    bed_dist_expansion = 0
  ) %>%
  dplyr::rename(ssr = bed_locus, locus = pos_locus)

missing_loci <-
  ssr_match_df_unneat %>%
  distinct(locus, .keep_all = T) %>%
  anti_join(select(ssr_match_df_by_start, locus))

ssr_match_df <-
  subset(ssr_match_df_unneat, locus %in% missing_loci$locus) %>%
  bind_rows(ssr_match_df_by_start)

yeast_genome <- read.fasta(genome_file)

ssr_bed_df_with_key <-
  ssr_bed_df_full %>%
  dplyr::rename(ssr = locus) %>%
  merge(ssr_match_df) %>%
  distinct() %>%
  dplyr::rename(SSR_start = start,
                SSR_end = end) %>%
  rowwise() %>%
  mutate(
    SSR_seq = toupper(paste(yeast_genome[[CHROM]][SSR_start:SSR_end],collapse = ''))
  ) %>%
  ungroup() %>%
  mutate(
    REPEAT_LEN = MOTIF_LEN*NUM_COPIES,
    AT_prop =
      str_count(SSR_seq,'A|T')/
      str_count(SSR_seq,'A|T|C|G'),
    AT_prop_approx = round(AT_prop*(AT_prop_approx_bins-1))/(AT_prop_approx_bins-1),
    short_ssr = REPEAT_LEN < long_ssr_min_length
  ) %>%
  mutate(`% A/T` = factor(
    case_when(AT_prop_approx == 0 ~ '<25%',
              AT_prop_approx == 0.5 ~ '25%-74%',
              AT_prop_approx == 1 ~ '>74%'),
    levels = c('>74%','25%-74%','<25%'))
  ) %>%
  mutate(MOTIF_NAME = factor(
    case_when(MOTIF_LEN==1 ~ 'homopolymer',
              MOTIF_LEN==2 ~ 'dinucleotide',
              MOTIF_LEN==3 ~ 'trinucleotide',
              MOTIF_LEN==4 ~ 'tetranucleotide',
              MOTIF_LEN > 4 ~ paste0(MOTIF_LEN,'-nucleotide')
              ),
    levels = c('homopolymer','dinucleotide','trinucleotide','tetranucleotide',paste0(5:max(MOTIF_LEN),'-nucleotide'))
    )
  ) %>%
  dplyr::select(-SSR_seq)

# identify mutations in repeats, telomeres, and centromeres
excluded_region_match_df <-
  recode_loci_overlap(
    pos_df,
    excluded_regions_df,
    alt_chrom_key = list('Mito' = 17*10^7),
    bed_dist_expansion = tel_cen_ltr_exclusion_dist
  )

################################################################################
# Process data
################################################################################
call_df <-
  curr_data %>%
  left_join(ssr_bed_df_with_key) %>%
  mutate(tel_cen_ltr = locus %in% excluded_region_match_df$pos_locus,
         call_allele_prop = AD/DP)

if (!is.na(additional_filter_exp)) {
  if (additional_filter_exp != 'NA'){
    call_df <-
      call_df %>%
      filter(!!rlang::parse_expr(additional_filter_exp))
  }
}

if (ploidy > 1){
  call_df <-
    call_df %>%
    mutate(
      allele = ifelse(grepl('\\.', allele), NA, allele),
      allele_1 = gsub('/.*','',allele),
      allele_2 = gsub('[A-z.,]+/','',allele),
      homozygous = allele_1 == allele_2
    ) %>%
    rename(GL_diff_ori = GL_diff) %>%
    mutate(GL_diff = ifelse(
      homozygous,
      GL_diff_ori/-log10(0.5),
      GL_diff_ori/-(0.5*(log10(p_qual)-log10(1-p_qual))-log10(0.5))
    ))
}

rel_call_df <- rel_call_fun(call_df, strain_anc_df)

################################################################################
# Identify non-motif mutations
################################################################################

# non-motif mutations are either SNMs (could be within SSR) or indels outside of
# SSR
non_motif_muts <-
  rel_call_df %>%
  filter(mutation) %>%
  filter(mut_bp_diff == 0 | is.na(ssr)) %>%
  mutate(superlocus = locus)

# loop, compressing superloci together (really probably only needs one loop)
superlocus_num = Inf
while (TRUE){
  # get list of mutated (super)loci
  non_motif_mut_candidates <-
    non_motif_muts %>%
    dplyr::select(superlocus, CHROM, start, end) %>%
    distinct() %>%
    group_by(superlocus, CHROM) %>%
    summarize(
      start = median(start),
      end = median(end)
    ) %>%
    ungroup() %>%
    dplyr::rename(locus = superlocus) %>%
    mutate(locus = as.factor(locus)) %>%
    arrange(CHROM, start)
  # find self-matches between (super)loci
  non_motif_self_match_df <-
    recode_loci_overlap(
      non_motif_mut_candidates,
      non_motif_mut_candidates,
      bed_dist_expansion = non_SSR_superlocus_dist
    ) %>%
    filter(pos_locus != bed_locus) %>%
    mutate(across(where(is.factor), as.character))
  if (length(unique(non_motif_self_match_df$pos_locus)) == 0){
    break
  } else{
    superlocus_num <- length(unique(non_motif_self_match_df$pos_locus))
    for (i in 1:nrow(non_motif_self_match_df)){
      curr_locus <- non_motif_self_match_df[i,'pos_locus']
      curr_superlocus <- non_motif_self_match_df[i,'bed_locus']
      non_motif_muts[
        non_motif_muts$superlocus == curr_locus,
        'superlocus'
        ] <- curr_superlocus
    }
  }
}

################################################################################
# Call SSR mutations
################################################################################
call_df_ssr <-
  call_df %>%
  filter(!tel_cen_ltr & !is.na(ssr))

GL_diff_thresh_df_slide_unfilt <- data.frame()
for (curr_window_size in slide_window_sizes){
  curr_GL_diff_thresh_df <-
    get_window_GL_diff_thresh(
      call_df_ssr,
      curr_window_size,
      quant_list
    ) %>%
    mutate(window_size = curr_window_size)
  GL_diff_thresh_df_slide_unfilt <-
    GL_diff_thresh_df_slide_unfilt %>%
    bind_rows(curr_GL_diff_thresh_df)
}

GL_diff_thresh_df_slide <-
  GL_diff_thresh_df_slide_unfilt %>%
  filter(ssr_locus_pop >= min_ssr_locus_per_category)


quant_thresh_genotypes_slide_cat <- data.frame()
quant_thresh_locus_counts <- data.frame()
for (w in unique(GL_diff_thresh_df_slide$window_size)){
  for (q in quant_list){
    curr_GL_diff_thresh_df_slide <-
      GL_diff_thresh_df_slide %>%
      filter(window_size == w & GL_diff_quant == q)
    curr_call_df <-
      call_df_ssr %>%
      inner_join(curr_GL_diff_thresh_df_slide) %>%
      filter(GL_diff < GL_diff_thresh)
    curr_mut_df <- rel_call_fun(curr_call_df, strain_anc_df)
    curr_quant_thresh_genotypes_slide_cat <-
      curr_mut_df %>%
      group_by_at(all_of(vars(ssr_grouping_cat))) %>%
      summarize(total_calls = sum(!is.na(mutation)),
                total_muts = sum(mutation, na.rm = T),
                ssr_loci = length(unique(ssr)),
                min_GL_diff = max(GL_diff)) %>%
      ungroup()  %>%
      mutate(window_size = w,
             GL_diff_quant = q)
    quant_thresh_genotypes_slide_cat <-
      quant_thresh_genotypes_slide_cat  %>%
      bind_rows(curr_quant_thresh_genotypes_slide_cat)
    curr_locus_mut_counts <-
      curr_mut_df %>%
      filter(mutation) %>%
      dplyr::count(locus) %>%
      mutate(window_size = w,
             GL_diff_quant = q)
    quant_thresh_locus_counts <-
      quant_thresh_locus_counts %>%
      bind_rows(curr_locus_mut_counts)
  }
}

################################################################################
# If call_df_ssr contains GL_diff_ori, use that to call mutations as well
################################################################################
ori_GL_ssr_mut_calls <- list()
if ('GL_diff_ori' %in% names(call_df_ssr)){
  call_df_ssr_ori <-
    call_df_ssr %>%
    select(-GL_diff) %>%
    rename(GL_diff = GL_diff_ori)
  
  GL_diff_thresh_df_slide_unfilt_ori <- data.frame()
  for (curr_window_size in slide_window_sizes){
    curr_GL_diff_thresh_df <-
      get_window_GL_diff_thresh(
        call_df_ssr_ori,
        curr_window_size,
        quant_list
      ) %>%
      mutate(window_size = curr_window_size)
    GL_diff_thresh_df_slide_unfilt_ori <-
      GL_diff_thresh_df_slide_unfilt_ori %>%
      bind_rows(curr_GL_diff_thresh_df)
  }
  
  GL_diff_thresh_df_slide_ori <-
    GL_diff_thresh_df_slide_unfilt_ori %>%
    filter(ssr_locus_pop >= min_ssr_locus_per_category)
  
  
  quant_thresh_genotypes_slide_cat_ori <- data.frame()
  quant_thresh_locus_counts_ori <- data.frame()
  for (w in unique(GL_diff_thresh_df_slide_ori$window_size)){
    for (q in quant_list){
      curr_GL_diff_thresh_df_slide <-
        GL_diff_thresh_df_slide_ori %>%
        filter(window_size == w & GL_diff_quant == q)
      curr_call_df <-
        call_df_ssr_ori %>%
        inner_join(curr_GL_diff_thresh_df_slide) %>%
        filter(GL_diff < GL_diff_thresh)
      curr_mut_df <- rel_call_fun(curr_call_df, strain_anc_df)
      curr_quant_thresh_genotypes_slide_cat <-
        curr_mut_df %>%
        group_by_at(all_of(vars(ssr_grouping_cat))) %>%
        summarize(total_calls = sum(!is.na(mutation)),
                  total_muts = sum(mutation, na.rm = T),
                  ssr_loci = length(unique(ssr)),
                  min_GL_diff = max(GL_diff)) %>%
        ungroup()  %>%
        mutate(window_size = w,
               GL_diff_quant = q)
      quant_thresh_genotypes_slide_cat_ori <-
        quant_thresh_genotypes_slide_cat_ori  %>%
        bind_rows(curr_quant_thresh_genotypes_slide_cat)
      curr_locus_mut_counts <-
        curr_mut_df %>%
        filter(mutation) %>%
        dplyr::count(locus) %>%
        mutate(window_size = w,
               GL_diff_quant = q)
      quant_thresh_locus_counts_ori <-
        quant_thresh_locus_counts_ori %>%
        bind_rows(curr_locus_mut_counts)
    }
  }
  ori_GL_ssr_mut_calls[['GL_diff_thresh_df_slide_unfilt']] <-
    GL_diff_thresh_df_slide_unfilt_ori
  ori_GL_ssr_mut_calls[['quant_thresh_genotypes_slide_cat']] <-
    quant_thresh_genotypes_slide_cat_ori
  ori_GL_ssr_mut_calls[['quant_thresh_locus_counts']] <- quant_thresh_locus_counts_ori
}


################################################################################
# Save output
################################################################################
save(
  call_df, rel_call_df, non_motif_muts, GL_diff_thresh_df_slide_unfilt,
  quant_thresh_genotypes_slide_cat, quant_thresh_locus_counts,
  ssr_bed_df_with_key, yeast_genome, ori_GL_ssr_mut_calls, excluded_regions_df,
  file = out_file
)

