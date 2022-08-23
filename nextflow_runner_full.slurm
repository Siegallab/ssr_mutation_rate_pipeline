#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=nf_full
#SBATCH --output=slurm_full.out

main_dir='/scratch/yp19/plavskin_seq_analysis'
nf_dir="${main_dir}/nf_scripts"

module purge
module load nextflow/20.10.0
cd $main_dir
nextflow run ${nf_dir}/sra_file_input_pipeline.nf -resume -with-timeline timeline_full.html -with-report report_full.html
