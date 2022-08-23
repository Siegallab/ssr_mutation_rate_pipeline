import sys
import math
import numpy as np
from Bio import pairwise2
import itertools
from copy import deepcopy


# Modified from analyze_overlaps.py by EP
# still use original 85% threshold for merging, but doesn't merge if motif length is same
# for overlapping regions that fail merging, split consecutive motifs down the middle of the overlap

def merge_regions(regions, merge_region_bool, location_array, nuc_array, min_start, match_score, mismatch_penalty, indel_penalty):
    merge_region_idx_list = np.where(merge_region_bool)[0]
    curr_seq_idx_list = np.any(location_array[merge_region_bool],axis = 0)
    curr_nuc_array = nuc_array[merge_region_bool][:,curr_seq_idx_list]
    curr_seq = ''.join(
        curr_nuc_array[(
            np.argmax(np.isin(curr_nuc_array,['A','T','C','G']),axis=0),
            range(0,curr_nuc_array.shape[1])
            )]
        )
    curr_motif = regions[merge_region_idx_list[0]][6]
    # calculate score using TRF parameters
    # this gives a score very close to TRF score
    curr_score = pairwise2.align.localms(
        curr_seq,''.join([curr_motif]*len(curr_seq)),
        match_score, mismatch_penalty, indel_penalty, indel_penalty,
        score_only = True)
    new_region = (
        regions[merge_region_idx_list[0]][0],
        min_start+np.min(np.where(curr_seq_idx_list)[0]),
        min_start+np.max(np.where(curr_seq_idx_list)[0]),
        len(curr_motif),
        float(len(curr_seq))/len(curr_motif),
        curr_score,
        curr_motif,
        curr_seq
        )
    regions = [r for i,r in enumerate(regions) if not (i in merge_region_idx_list)]
    regions.append(new_region)
    return(regions)

def get_nuc_array(regions, region_size, min_start):
    nuc_array=np.chararray((len(regions),region_size),unicode=True)
    for index,region in enumerate(regions):
        curr_start = region[1]
        curr_stop = region[2]
        curr_nucs = np.array(list(region[7]))
        nuc_array[index,(curr_start-min_start):(curr_stop-min_start+1)] = curr_nucs
    return(nuc_array)

def get_possible_stretches(start, region_end, min_ssr_len):
    # calculates possible stretches of the given region from start to 
    # region_end that result in a length of either 0 or min_ssr_len
    poss_ends = [start]+list(range(start+min_ssr_len-1,region_end+1))
    stretch_tuples = list(zip([start]*len(poss_ends),poss_ends))
    return(stretch_tuples)

def get_region_start_end_tuples(region, min_ssr_len):
    poss_starts = range(region[1],region[2]+1)
    start_end_tuples = list(itertools.chain(*[get_possible_stretches(s, region[2], min_ssr_len) for s in poss_starts]))
    return(start_end_tuples)

def get_min_ssr_len(region, min_homopolymer_rep, min_nonhomopolymer_rep):
    curr_motif_len = region[3]
    if curr_motif_len == 1:
        min_ssr_len = min_homopolymer_rep
    else:
        min_ssr_len = int(curr_motif_len*min_nonhomopolymer_rep)
    return(min_ssr_len)


def get_optimal_score_split(curr_regions, min_homopolymer_rep, min_nonhomopolymer_rep, min_start, max_stop, min_score_per_bp, match_score, mismatch_penalty, indel_penalty, include_length_constraint = False):
    if len(curr_regions) == 2:
        # try removing one region of the other, or breaking them up in 
        # such a way that all breakpoints are at shared positions
        # this may lost sequences whose overall match to motif is very 
        # poor but whose subsequence has a strong match to motif
        # (but speeds up code a lot)
        nuc_array = get_nuc_array(curr_regions, (max_stop+1-min_start), min_start)
        location_array = np.isin(nuc_array,['A','T','C','G'])
        intersect_bool = np.sum(location_array,axis = 0) == 2
        first_intersect_idx = np.where(intersect_bool)[0][0]
        last_intersect_idx = np.where(intersect_bool)[0][-1]
        start_end_tuple_list = []
        for i,region in enumerate(curr_regions):
            region_start_idx = np.where(location_array[i])[0][0]
            region_end_idx = np.where(location_array[i])[0][-1]
            start_end_tuples = []
            min_ssr_len = get_min_ssr_len(region, min_homopolymer_rep, min_nonhomopolymer_rep)
            if (region_start_idx < first_intersect_idx):
                # start at beginning of region, ends in or right before intersect
                poss_starts = region_start_idx+min_start
                poss_ends_unfilt = [first_intersect_idx-1+min_start]+list(np.where(intersect_bool)[0]+min_start)
                poss_ends = [poss_starts]+[e for e in poss_ends_unfilt if ((e-poss_starts) >= (min_ssr_len-1))]
                start_end_tuples = start_end_tuples+list(zip([poss_starts]*len(poss_ends),poss_ends))
            if (region_end_idx > last_intersect_idx):
                # start inside intersect, end at end of region
                poss_starts_unfilt = range(first_intersect_idx+min_start, last_intersect_idx+min_start+2)
                poss_ends = region_end_idx+min_start
                poss_starts = [poss_ends]+[s for s in poss_starts_unfilt if ((poss_ends-s) >= (min_ssr_len-1))]
                start_end_tuples = start_end_tuples+list(zip(poss_starts,[poss_ends]*len(poss_starts)))
            if ((region_start_idx < first_intersect_idx) & (region_end_idx > last_intersect_idx)):
                # overlap is contained entirely within region, include full-region option
                start_end_tuples = start_end_tuples+[tuple((region_start_idx+min_start, region_end_idx+min_start))]
            if ((region_start_idx == first_intersect_idx) & (region_end_idx == last_intersect_idx)):
                # region entirely overlapped, anything can be start/end
                start_end_tuples = get_region_start_end_tuples(region, min_ssr_len)
            start_end_tuple_list.append(start_end_tuples)
    else:
        # try all possible start and stop values for every region
        # find combo that maximizes score (and, among max scores, 
        # minimizes region number)
        start_end_tuple_list = []
        for region in curr_regions:
            min_ssr_len = get_min_ssr_len(region, min_homopolymer_rep, min_nonhomopolymer_rep)
            start_end_tuples = get_region_start_end_tuples(region, min_ssr_len)
            start_end_tuple_list.append(start_end_tuples)
    start_end_tuple_combos_unfilt = list(itertools.product(*start_end_tuple_list))
    # only calculate scores if curr_regions non-overlapping and non-blank
    # also don't calculate scores for curr_regions that comprise < 80% 
    # of original region length; unlikely to have higher scores
    # (for speed)
    start_end_tuple_combos = []
    if include_length_constraint:
        min_region_length = 0.5*(max_stop+1-min_start)
    else:
        min_region_length = 0
    for tuple_combo in start_end_tuple_combos_unfilt:
        tuple_combo_point_list = list(itertools.chain(*[list(range(s,t+1)) for s,t in tuple_combo]))
        # only add if length is > min_region_length and no position covered twice (no overlap)
        if not (
            len(tuple_combo_point_list) < min_region_length or
            len(tuple_combo_point_list) != len(set(tuple_combo_point_list))):
            start_end_tuple_combos.append(tuple_combo)
    max_score_counter = 0
    min_len = len(curr_regions)
    new_regions = curr_regions
    nuc_array = get_nuc_array(curr_regions, (max_stop+1-min_start), min_start)
    for j,start_end_tuple in enumerate(start_end_tuple_combos):
        curr_score = 0
        curr_len = 0
        curr_new_regions = []
        for i,region in enumerate(curr_regions):
            r_start = start_end_tuple[i][0]
            r_stop = start_end_tuple[i][1]
            curr_seq = ''.join(nuc_array[i,(r_start-min_start):(r_stop+1-min_start)])
            min_ssr_len = get_min_ssr_len(region, min_homopolymer_rep, min_nonhomopolymer_rep)
            if len(curr_seq) >= min_ssr_len:
                curr_len = curr_len+1
                curr_perfect_seq = ''.join(region[6]*int(math.ceil(len(curr_seq)/region[3])))
                subscore = \
                    pairwise2.align.localms(
                        curr_seq,curr_perfect_seq,
                        match_score, mismatch_penalty, indel_penalty, indel_penalty,
                        score_only = True)
                if (subscore == []):
                    subscore = 0
                subscore_per_bp = float(subscore)/len(curr_seq)
                if subscore_per_bp >= min_score_per_bp:
                    curr_score = curr_score + subscore
                    curr_new_regions.append((
                        region[0],
                        r_start,
                        r_stop,
                        region[3],
                        float(len(curr_seq))/region[3],
                        subscore,
                        region[6],
                        curr_seq
                        ))
        if (curr_score > max_score_counter) or ((curr_score == max_score_counter) and (curr_len < min_len)):
            new_regions = list(curr_new_regions)
            max_score_counter = curr_score
            min_len = curr_len
    return(new_regions)

def main():
    data       = open(sys.argv[1], "r")
    merge_pass = open(sys.argv[2], "w")
    merge_fail = open(sys.argv[3], "w")
    min_homopolymer_rep = int(sys.argv[4])
    min_nonhomopolymer_rep = int(sys.argv[5])
    match_score = int(sys.argv[6])
    mismatch_penalty = -int(sys.argv[7])
    indel_penalty = -int(sys.argv[8])
    min_score_per_bp = float(sys.argv[9])
    pass_count = 0
    fail_count = 0
    regions    = []
    cur_chrom  = ""
    min_start  = 0
    max_stop   = -1 

    for line in data:
        tokens              = line.strip().split("\t")
        chrom               = tokens[0]
        start, stop, period = map(int, tokens[1:4])
        motif               = tokens[4]
        num_repeats         = float(tokens[5])
        score               = float(tokens[6])
        sequence            = tokens[7]
        # here, we're looping through SORTED repeats until we hit a 
        # different chromosome or until the start of the next repeat 
        # is past the stop of the most recent ones (i.e. 
        # non-overlapping)
        if chrom != cur_chrom or start > max_stop:
            if len(regions) > 0:
                # Analyze the previous set of overlapping regions
                region_size = max_stop - min_start + 1
                max_score   = 0
                max_index   = -1
                cov_frac    = 0.0
                # merge any overlapping regions with identical motifs
                # there may be some clever way to do this all at once, but 
                # instead looping through and merging best merge candidates
                if len(regions) > 1:
                    motif_merges_remaining = True
                    # check that region length is equal to length of sequence for each region, merge/split code below fails
                    regions = [r for r in regions if (len(r[7]) == (r[2]-r[1]+1))]
                while motif_merges_remaining:
                    nuc_array = get_nuc_array(regions, region_size, min_start)
                    location_array = np.isin(nuc_array,['A','T','C','G'])
                    curr_overlap_mat = np.matmul(location_array,np.transpose(location_array))
                    motif_list = [r[6] for r in regions]
                    motif_identity_mat = np.array(motif_list) == np.transpose(np.array(motif_list)[np.newaxis])
                    motif_merge_array = np.logical_and(motif_identity_mat, curr_overlap_mat)
                    motif_merge_scores = np.sum(motif_merge_array.astype(int),axis = 1)
                    if np.any(motif_merge_scores > 1):
                        merge_region_bool = motif_merge_array[np.argmax(motif_merge_scores)]
                        regions = merge_regions(regions, merge_region_bool, location_array, nuc_array, min_start, match_score, mismatch_penalty, indel_penalty)
                    else:
                        motif_merges_remaining = False
                # merge any overlaps where one array is a subarray of another *and* has a lower score (and lower score per bp)
                # (likely unnecessary step with new get_optimal_score_split)
                if len(regions) > 1:
                    full_overlap_merges_remaining = True
                else:
                    full_overlap_merges_remaining = False
                while full_overlap_merges_remaining:
                    nuc_array = get_nuc_array(regions, region_size, min_start)
                    location_array = np.isin(nuc_array,['A','T','C','G'])
                    curr_overlap_mat = np.matmul(location_array,np.transpose(location_array))
                    # find sequences that overlap only a single sequence *and* whose 
                    # overlap is with that sequence for its entire span (i.e. it always
                    # overlaps it)
                    single_overlap_bool = np.sum(curr_overlap_mat,axis = 0)==2
                    fully_contained_bool = \
                        np.all(np.isin((np.sum(location_array,axis=0)*location_array),[0,2]),axis=1)
                    fully_contained_single_overlap_bool = \
                        np.logical_and(single_overlap_bool,fully_contained_bool)
                    curr_overlap_mat_nonself = curr_overlap_mat.copy()
                    np.fill_diagonal(curr_overlap_mat_nonself,0)
                    # identify index of parent region for each region in fully_contained_single_overlap_bool
                    parent_region_idx = \
                        np.where(curr_overlap_mat_nonself*fully_contained_single_overlap_bool)
                    # initialize parent_array with non-real region index
                    parent_array = np.array([len(regions)+1]*len(regions))
                    parent_array[parent_region_idx[1]] = parent_region_idx[0]
                    # identify regions with lower score and per-bp score than the bigger region they're part of
                    score_array = np.array([r[5] for r in regions])
                    score_per_bp_array = np.array([float(r[5])/len(r[7]) for r in regions])
                    parent_score_array = np.zeros_like(score_array)
                    parent_score_array[parent_region_idx[1]] = score_array[parent_region_idx[0]]
                    parent_score_per_bp_array = np.zeros_like(score_array)
                    parent_score_per_bp_array[parent_region_idx[1]] = score_per_bp_array[parent_region_idx[0]]
                    subregion_merge_bool = np.logical_and(
                        parent_score_array>score_array,
                        parent_score_per_bp_array>score_per_bp_array,
                        fully_contained_single_overlap_bool
                        )
                    # CAN ONLY MERGE WITH ONE PARENT AT A TIME WITHOUT BREAKING INDEXING
                    merge_parent_candidates = np.unique(parent_array[subregion_merge_bool])
                    if len(merge_parent_candidates)==0:
                        full_overlap_merges_remaining = False
                    else:
                        curr_merge_parent = merge_parent_candidates[0]
                        curr_region_parent_array = parent_array==curr_merge_parent
                        merge_region_bool = np.logical_and(
                            subregion_merge_bool, curr_region_parent_array
                            )
                        merge_region_bool[curr_merge_parent] = True
                        regions = merge_regions(regions, merge_region_bool, location_array, nuc_array, min_start, match_score, mismatch_penalty, indel_penalty)
                if any([r[5]>0 for r in regions]):
                    if len(regions) == 1:
                        region = regions[0]
                        merge_pass.write("%s\t%s\t%d\t%d\t%.1f\t%s\n"%(region[0], region[1], region[2], region[3], region[4], region[6]))
                        pass_count += (len(regions) > 1)
                    else:
                        # find split of regions that optimizes score
                        # if more than two regions, too many possibilities
                        # optimize single overlap at a time
                        opt_merges_remain = True
                        while opt_merges_remain:
                            region_scores = np.array([r[5] for r in regions])
                            region_score_sum_mat = region_scores + np.transpose(region_scores[np.newaxis])
                            nuc_array = get_nuc_array(regions, region_size, min_start)
                            location_array = np.isin(nuc_array,['A','T','C','G'])
                            curr_overlap_mat = np.matmul(location_array,np.transpose(location_array))
                            curr_overlap_mat_nonself = curr_overlap_mat.copy()
                            np.fill_diagonal(curr_overlap_mat_nonself,0)
                            overlap_nums = np.sum(curr_overlap_mat_nonself, axis = 1)
                            if np.any(overlap_nums>0):
                                # optimize merge of two overlapping regions with biggest combined score first
                                masked_region_score_sum_mat = np.ma.masked_array(
                                    region_score_sum_mat, mask = curr_overlap_mat_nonself==0
                                    )
                                region_idx_tuple = np.unravel_index(
                                    np.argmax(masked_region_score_sum_mat),
                                    masked_region_score_sum_mat.shape
                                    )
                                curr_regions_to_split = [r for i,r in enumerate(regions) if i in region_idx_tuple]
                                curr_unsplit_regions = [r for i,r in enumerate(regions) if i not in region_idx_tuple]
                                curr_split_regions = get_optimal_score_split(curr_regions_to_split, min_homopolymer_rep, min_nonhomopolymer_rep, min_start, max_stop, min_score_per_bp, match_score, mismatch_penalty, indel_penalty)
                                regions = curr_split_regions + curr_unsplit_regions
                            else:
                                opt_merges_remain = False
                        for region in regions:
                            merge_fail.write("%s\t%s\t%d\t%d\t%.1f\t%s\n"%(region[0], region[1], region[2], region[3], region[4], region[6]))
                            fail_count += (len(regions) > 1)
            # Reset the region info
            regions   = []
            min_start = start
            max_stop  = stop
            cur_chrom = chrom
        else:
            max_stop = max(max_stop, stop)
        regions.append((chrom, start, stop, period, num_repeats, score, motif, sequence))

    merge_pass.close()
    merge_fail.close()
    data.close()

if __name__ == "__main__":
    main()
