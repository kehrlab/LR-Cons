library(tidyverse)

# shared variables and paths
output_dir <- "results/plots_final/"
numbers_dir <- "results/numbers_final/"

# specific to real data
variant_file_real <- "data/references/HG002_assembly/variant_regions_500.txt"
results_dir_real <- "results/real/"
annotation_real <- "data/references/t2t-chm13/variant_annotation.tsv"
stats_ont_real <- "data/real/HG002/variant_region_stats_ont_32.51.tsv"
stats_pacbio_real <- "data/real/HG002/variant_region_stats_pacbio_32.51.tsv"
variant_coverage_whitelist <- "results/cov_var_list.tsv"
variant_coverage_whitelist_old <- "results/cov_var_list_with_old.tsv"
variant_overlap_whitelist_500 <- "data/references/HG002_assembly/non_overlap_variants_500.txt"
variant_overlap_whitelist_1000 <- "data/references/HG002_assembly/non_overlap_variants_1000.txt"

# specific to GIAB simulation
variant_file_giab <- "data/simulation/giab/variant_regions_500_0.tsv"
results_dir_giab <- "results/giab_simulation/"
giab_overlap_whitelist <- "data/simulation/giab/non_overlap_variants_500.txt"


# specific to uncertainty simulation
results_dir_sim <- "results/uncertainty_simulation/"

# specific to main simulation
variant_file <- "data/simulation/variant_regions_500_0.tsv"
results_dir <- "results/simulation"

# specific to kmer-size simulation
variant_file_kmeans <- "data/simulation/variant_regions_500_0.tsv"
results_dir_kmeans <- "results/kmeans_simulation/"


# for plotting
my_colors <- c(
  '#ffd266ff', # gelb
  '#f58080ff', #rot
  '#6a00e5ff',
  '#12beccff',
  '#fcab79ff',
  '#0059ecff',
  '#000c52',
  "lightgreen", # not LIT
  '#ff69b4' # not LIT
)

method_order <- c(
  "abPOA",
  "Multi POA",
  "RapidFuzz",
  "k-mer clust",
  "A priori",
  "True phase"
)

method_order_v2 <- c(
  "abPOA",
  "Multi POA",
  "RapidFuzz",
  "k-mer clust",
  "A priori",
  "True phase",
  "Multi POA v2"
)

method_order_real <- c(
  "abPOA",
  "Multi POA",
  "RapidFuzz",
  "k-mer clust",
  "A priori"
)


group_order <- c(4, 5, 2, 7, 3, 6, 1)
group_names <- c(
  "HOM_ONE", "HOM_TWO", "HOM_NONE",
  "HET_DISTINCT", "HET_ONE", "HET_REDUNDANT", "HET_NONE"
)

custom_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white")
  )

custom_theme_poster <- theme_bw() +
  theme(
    plot.title = element_text(size = 26),
    plot.subtitle = element_text(size = 22),
    legend.title = element_text(size = 24), 
    legend.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 22),
    strip.background = element_rect(fill = "white")
  )