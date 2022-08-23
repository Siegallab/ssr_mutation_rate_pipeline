# Functions to identify overlapping loci

library(Matrix)

make_chrom_pos <- function(chrom_list, alt_chrom_key = NA, mult_val = 10^7){
  abs_pos_list <- suppressWarnings(as.numeric(as.roman(chrom_list))*mult_val)
  if (!all(is.na(alt_chrom_key))){
    for (key in names(alt_chrom_key)){
      abs_pos_list[chrom_list == key] <- alt_chrom_key[[key]]
    }
  }
  return(abs_pos_list)
}

create_pos_mat <- function(mat_df, col_num){
  # make a sparse matrix based on start and end positions in mat_df
  mat_df <-
    mat_df %>%
    mutate(width = end - start + 1) %>%
    rowid_to_column()
  # vectorized 'seq' from https://stackoverflow.com/a/15917618/8082611
  seq_vect <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  out_mat <- sparseMatrix(
    i = rep(mat_df$rowid, mat_df$width),
    j = unlist(seq_vect(mat_df$start, mat_df$end)),
    x = 1,
    dims = c(nrow(mat_df), col_num)
  )
  return(out_mat)
}

recode_loci_overlap <- function(
  position_df, bed_df, bed_dist_expansion = 0, alt_chrom_key = list()
){
  # identify loci in bed_df that overlap loci in pos_df
  # if bed_dist_expansion > 0, it's added to start and end of each locus in bed_df
  # add 'absolute position' (chromosome-independent) to dataframes
  pos_mat_df <-
    position_df %>%
    mutate(
      start = make_chrom_pos(CHROM, alt_chrom_key = alt_chrom_key)+start,
      end = make_chrom_pos(CHROM, alt_chrom_key = alt_chrom_key)+end
    ) %>%
    select(c(start, end))
  bed_mat_df <-
    bed_df %>%
    mutate(
      mod_start = 
        make_chrom_pos(CHROM, alt_chrom_key = alt_chrom_key) +
        start -
        bed_dist_expansion,
      end = 
        make_chrom_pos(CHROM, alt_chrom_key = alt_chrom_key) +
        end +
        bed_dist_expansion,
      start = 
        ifelse(mod_start < 0, 0, mod_start)
    ) %>%
    select(c(start, end))
  # create sparse matrix of positions where SSRs are
  max_genome_pos = max(c(bed_mat_df$end, pos_mat_df$end))
  bed_mat <- create_pos_mat(bed_mat_df, max_genome_pos)
  pos_mat <- create_pos_mat(pos_mat_df, max_genome_pos)
  # make overlap matrix with columns (j) corresponding to bed_df rows,
  # rows (i) corresponding to pos_df rows, and 
  # the number in the position (x) corresponding
  # to the number of overlapping positions between the two
  overlap_mat <- as(pos_mat %*% t(bed_mat),"dgTMatrix")
  # sparse matrices are 0-based apparently?
  # see here https://slowkow.com/notes/sparse-matrix/
  # this indexing only returns non-zero positions
  match_df_locus <- data.frame(
    pos_locus = position_df$locus[overlap_mat@i+1],
    bed_locus = bed_df$locus[overlap_mat@j+1]
  )
  pos_multimatch_df <-
    match_df_locus %>%
    dplyr::count(pos_locus) %>%
    filter(n > 1)
  bed_multimatch_df <-
    match_df_locus %>%
    dplyr::count(bed_locus) %>%
    filter(n > 1)
  pos_multimatch_num <- nrow(pos_multimatch_df)
  bed_multimatch_num <- nrow(bed_multimatch_df)
  if (pos_multimatch_num > 0){
    print(paste(c("Warning, multiple matches in",
                  pos_multimatch_num,
                  "of the following loci:",
                  as.character(pos_multimatch_df$pos_locus[1:min(pos_multimatch_num,100)])
    ), collapse = " "))
  }
  if (bed_multimatch_num > 0){
    print(paste(c("Warning, multiple matches in",
                  bed_multimatch_num,
                  "of the following loci:",
                  as.character(bed_multimatch_df$bed_locus[1:min(bed_multimatch_num,100)])
    ), collapse = " "))
  }
  return(match_df_locus)
}

# Functions to call relative mutations

haploid_ma_rel_call_fun <- function(call_df, strain_anc_df){
  
  parent_call_df <-
    call_df %>%
    select(
      c(
        strain, locus, GT, GL_diff, DP, allele, call_allele_prop
      )
    ) %>%
    filter(strain %in% strain_anc_df$anc_strain) %>%
    rename_with(.cols=-c(locus), function(x){paste0("anc_",x)}) %>%
    full_join(strain_anc_df)
  
  rel_call_df <-
    call_df %>%
    filter(!(strain %in% strain_anc_df$anc_strain)) %>%
    inner_join(parent_call_df) %>%
    mutate(mutation = GT != anc_GT,
           min_DP = pmin(DP, anc_DP),
           min_GL_diff = pmax(GL_diff, anc_GL_diff)
    ) %>%
    filter(!tel_cen_ltr) %>%
    mutate(
      mut_bp_diff = ifelse(
        allele != anc_allele,
        nchar(allele)-nchar(anc_allele),
        NA
      )
    )
  return(rel_call_df)
}

diploid_ma_rel_call_fun <- function(call_df, strain_anc_df){
  parent_call_df <-
    call_df %>%
    select(
      c(
        strain, locus, GT, GL_diff, DP, allele_1, allele_2, call_allele_prop
      )
    ) %>%
    filter(strain %in% strain_anc_df$anc_strain) %>%
    rename_with(.cols=-c(locus), function(x){paste0("anc_",x)}) %>%
    full_join(strain_anc_df)
  
  rel_call_df <-
    call_df %>%
    filter(!(strain %in% strain_anc_df$anc_strain)) %>%
    inner_join(parent_call_df) %>%
    mutate(mutation = GT != anc_GT,
           min_DP = pmin(DP, anc_DP),
           min_GL_diff = pmax(GL_diff, anc_GL_diff)
    ) %>%
    mutate(mutation = ifelse(grepl('\\.', GT), NA, mutation)) %>%
    filter(!tel_cen_ltr) %>%
    mutate(
      mut_bp_diff = case_when(
        (allele_1 != anc_allele_1) & (allele_2 == anc_allele_2) ~ nchar(allele_1)-nchar(anc_allele_1),
        (allele_1 == anc_allele_1) & (allele_2 != anc_allele_2) ~ nchar(allele_2)-nchar(anc_allele_2)
      )
    )
  return(rel_call_df)
}

# Functions for calling SSR mutations

get_window_GL_diff_thresh <- function(call_df, copy_num_window_size, quantiles){
  one_side_window_size <- copy_num_window_size/2
  copy_nums <- sort(unique(call_df$NUM_COPIES))
  out_df <- data.frame()
  for (n in copy_nums){
    window_start <- n - one_side_window_size
    window_end <- n + one_side_window_size
    out_df <-
      call_df %>%
      filter(NUM_COPIES > window_start & NUM_COPIES <= window_end) %>%
      group_by(MOTIF_LEN, AT_prop_approx, `% A/T`) %>%
      summarize(GL_diff_quant = list(quantiles),
                GL_diff_thresh = list(-quantile(-GL_diff, quantiles, names = F, na.rm = T)),
                pop = sum(!is.na(GL_diff)),
                ssr_locus_pop = length(unique(ssr))) %>%
      unnest(c(GL_diff_quant, GL_diff_thresh)) %>%
      mutate(NUM_COPIES = n) %>%
      bind_rows(out_df)
  }
  return(out_df)
}


