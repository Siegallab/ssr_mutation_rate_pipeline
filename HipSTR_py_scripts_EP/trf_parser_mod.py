import sys
import numpy as np

num_tokens     = 17
score_index    = 7
period_index   = 4  # Actually the size of the consensus sequence, but this is more accurate
motif_index    = 13
sequence_index = 14
start_index    = 0
stop_index     = 1
nrepeat_index  = 3

def min_perm(seq):
    min_perm = seq
    for i in xrange(len(seq)):
        other = seq[i:]+seq[0:i]
        if other < min_perm:
            min_perm = other
    return min_perm

def rev_complement(seq):
    rev     = seq[::-1]
    mapping = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    res     = ""
    for i in xrange(len(rev)):
        res = res+mapping[rev[i]]
    return res

def canonical_motif(motif):
    return min(min_perm(motif), min_perm(rev_complement(motif)))

def create_filtered_trf_bed_file(input_files, max_period, min_hom_rep, min_non_hom_rep):
    period_vals       = np.arange(1,(max_period+1))
    period_thresholds = np.append([min_hom_rep*2], period_vals[1:]*2*min_non_hom_rep)
    filt_count        = 0
    keep_count        = 0

    for i in range(len(input_files)):
        data  = open(input_files[i], "r")
        chrom = data.readline().strip()[1:]
        for line in data:
            tokens = line.strip().split()
            period = int(tokens[period_index])
            if period > max_period:
                filt_count += 1
                continue

            score = int(tokens[score_index])
            if score < period_thresholds[period-1]:
                filt_count += 1
                continue

            new_tokens    = [chrom] + list(map(lambda x: tokens[x], [start_index, stop_index, period_index, motif_index, nrepeat_index, score_index, sequence_index]))
            new_tokens[4] = min_perm(new_tokens[4])
            print("\t".join(new_tokens))
            keep_count += 1

        data.close()
    sys.stderr.write("Kept %d out of %d records (%.2f%%) with sufficiently high scores\n"%(keep_count, keep_count+filt_count, 100.0*keep_count/(keep_count+filt_count)))
    return

def main():
    files = sys.argv[1].strip().split(",")
    max_period = int(sys.argv[2])
    min_hom_rep = int(sys.argv[3])
    min_non_hom_rep = int(sys.argv[4])

    create_filtered_trf_bed_file(files, max_period, min_hom_rep, min_non_hom_rep)

if __name__ == "__main__":
    main()
