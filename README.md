# Repository for "Evaluating local assembly strategies for long-read structural variant valling"

Scripts, code and Snakemake workflow for simulation of SVs on T2T-CHM13 and ONT long reads
for evaluation of different ways to generate properly phased local consensus sequences
at SV loci.

## Project structure
Scripts for data generation (python) and evaluation (R) are located in `scripts/`. 

C++ source code for extraction of region sequences* from alignment data,
creation of consensus sequences with different approaches and alignment of consensus sequences
against true allele sequences for all variants is located in `src/`. 
Code can be compiled with `make all`. Dependencies are `SPOA` and `seqan2`.

A conda environment for compilation can be installed with
```
conda create -n lr_benchmark_env -c conda-forge -c bioconda htslib=1.21 seqan=2.4 spoa=4.1.5 cxx-compiler libdeflate=1.22 zlib rapidfuzz-cpp abpoa=1.5.6
```

## Results 
Plots underlying the Figures in the paper were generated with [`scripts/create_paper_plots.R`](./scripts/create_paper_plots.R) and can be found in [`results/plots`](./results/plots).
Supplementary Plots were generated with [`scripts/create_supplementary_figures.R`](scripts/create_supplementary_figures.R) and can be found in [`results/supplementary_plots`](./results/supplementary_plots).
Data shown in the plots or mentioned in the paper can be found as tables in [`results/numbers`](./results/numbers).

Relevant data for the plot generation can be found in:
  - `results/ont` and `results/pacbio`, divided in subfolders containig output of `./evaluate_consensus sequences`: tables with results of consensus sequence alignments to true haplotype sequences
  - `data/time_data`, containing tables with statistics about consensus sequence generation methods, including runtime for each locus

Simulated data and generated consensus sequences are not contained in this repository.

## Current Workflow
![image](dag.svg)

1. Insertion of 16000 SVs on chromosome 1 of T2T-CHM13.
    - Deletions
    - Inversions
    - Tandem Duplications
    - Compound Heterozygous Deletions
    - Compound Heterozygous Inv + Del (no overlap)
    Simulated sizes: 50, 200, 500, 1000
2. Simulation of long read data using PBSIM
    - ONT, PacBio HiFi
    - 20x, 30x, 60x coverage
    - Accuracy 99%
3. Alignment with Minimap2
4. Extraction of sequences spanning the variant regions from long read alignments
5. Creation of consensus sequence(s) using 6 different methods
    - Heaviest Bundle Algorithm (abPOA)
    - K-Mer based clustering + POA (SPOA)
    - Iterative generation of multiple POA graphs (SPOA)
    - A priori phasing
    - RapidFuzz 
    - Based on haplotype of origin (known through simulation)
6. Alignment of consensus sequences against true allele sequences for each variant
