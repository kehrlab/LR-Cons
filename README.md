# Repository for "Evaluating local assembly strategies for long-read structural variant valling in human genomes"

Scripts, code and Snakemake workflow 
for evaluation of different ways to generate properly phased local consensus sequences
at SV loci.
Includes download and preprocessing of real long read sequencing data, as well as simulation of SVs and long reads.

## Project structure
Scripts for data generation (python) are located in `scripts/` and scripts for evaluation (R) are located in `evaluation_scripts/`. 

C++ source code for extraction of region sequences* from alignment data,
creation of consensus sequences with different approaches and alignment of consensus sequences
against true allele sequences for all variants is located in `src/`. 
Code can be compiled with `make all`. Dependencies are `SPOA`, `seqan2`, `rapidfuzz` and `abPOA`.

A conda environment for compilation can be installed with
```
conda create -n lr_benchmark_env -c conda-forge -c bioconda htslib=1.21 seqan=2.4 spoa=4.1.5 cxx-compiler libdeflate=1.22 zlib rapidfuzz-cpp abpoa=1.5.6
```

## Results 
Plots underlying the Figures in the paper were generated with [`evaluation_scripts/create_plots.R`](./evaluation_scripts/create_plots.R) and can be found in [`results/plots_final`](./results/plots_final).
Data shown in the plots or mentioned in the paper can be found as tables in [`results/numbers_final`](./results/numbers_final).

Relevant data for the plot generation can be mainly found in `results` and is divided in subfolders containing output of `./evaluate_consensus sequences` (tables with results of consensus sequence alignments to true haplotype sequences):
  - `results/real`: generated with `snakemake real` and `snakemake real_old`
  - `results/simulation`: generated with `snakemake sim`
  - `results/giab_simulation`: generated with `snakemake sim_giab`
  - `results/kmeans_simulation`: generated with `snakemake benchmark_k`
  - `results/uncertainty_simulation`: generated with `snakemake uncertainty`

Runtimes and other metrics about created consensus sequences are contained in `data/time_benchmark` and can be generated with `snakemake time_benchmark`.

Simulated data and generated consensus sequences are not contained in this repository.
Real data will be automatically downloaded and stored in `data/`, unless otherwise specified in `config.yaml`.
