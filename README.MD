# Substitution and SSR mutation calling pipeline

This pipeline represents the workflow used to call mutations in wild type and *Msh3*-deficient *Saccharomyces cerevisiae* strains in the Plavskin, de Biase, *et al.* preprint.

The basic structure of the pipeline is described in the preprint.

In this pipeline, the sequencing data is pulled from the Sequence Read Archive.

## Requirements

Required software is listed in the [config file](nextflow.config). Note that the pipeline assumes module is installed. The python scripts require biopython, pyvcf, and numpy to be installed in the local environment.

The pipeline is currently configured to be run through slurm, but this can be changed in the [config file](nextflow.config).

## Running the pipeline

The pipeline can be submitted as a batch job on slurm using the [submission script](nextflow_runner_full.slurm). It runs [sra_file_input_pipeline.nf](nf_scripts/sra_file_input_pipeline.nf), which pulls *fastq* files from SRA and launches [generic_workflow.nf](nf_scripts/generic_workflow.nf). This in turn creates the key output files, including a .csv file with genotyping data for every position in the genome with a possible substitution, and every SSR (regardless of the probability that the latter contains an indel).

After the pipeline has been run, figures and key results from the preprint can be reproduced with the [SSR_mutation_rate_results_plotting.Rmd R notebook](SSR_mutation_rate_results_plotting.Rmd).