import itertools
import sys
import pandas as pd

fasta_output = sys.argv[1]

max_length = 10

def expand_grid(data_dict):
	rows = itertools.product(*data_dict.values())
	return pd.DataFrame.from_records(rows, columns=data_dict.keys())

nucleotides = ['A','T','C','G']
motif_df = pd.DataFrame()
for per_length in [1,2,3]:
	nuc_combo = itertools.product(nucleotides, repeat=per_length)
	if per_length == 1:
		curr_list = list(''.join(i) for i in nuc_combo)
	else:
		curr_list = list(''.join(i) for i in nuc_combo if len(set(i)) != 1)
	curr_df = expand_grid({'motif':curr_list, 'motif_len':[per_length], 'seq_length':range(4,max_length+1)})
	motif_df = pd.concat([motif_df, curr_df])

motif_df = (motif_df.
	query('seq_length >= motif_len * 2').
	assign(
		full_seq=lambda x: x.motif*max_length+x.motif.str[0],
		num_copies=lambda x: (x.seq_length/x.motif_len).round(2),
		exclude_before=lambda x: x.motif.str[-1]
		).
	assign(
		ssr_seq=lambda x: x.apply(
			lambda row: row['full_seq'][0:row['seq_length']], axis = 1
			),
		exclude_after=lambda x: x.apply(
			lambda row: row['full_seq'][row['seq_length']], axis = 1
			)
		).
	assign(
		fasta_name = lambda x: x.apply(
			lambda row: '>'+'-'.join([row['motif'],str(row['motif_len']),str(row['num_copies']),str(2*row['seq_length']),row['ssr_seq']]), axis = 1
			),
		fasta_pattern = lambda x: '[^'+x.exclude_before+']('+x.ssr_seq+')[^'+x.exclude_after+']'
		)
	)

with open(fasta_output, "w") as fasta_file:
	fasta_file.write('\n'.join(list(itertools.chain(*zip(motif_df.fasta_name,motif_df.fasta_pattern)))))

