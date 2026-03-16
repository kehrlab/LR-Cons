library(tidyverse)
source("scripts/misc.R")

# -- function definitions --- #
determine_group <- function(HOMOZYGOUS, N, N_unique, NoSeq)
{
  groups <- rep(NA, length(HOMOZYGOUS))
  
  groups[!is.na(HOMOZYGOUS) & N == 1 & NoSeq & !HOMOZYGOUS] <- "G"
  groups[!is.na(HOMOZYGOUS) & N == 1 & NoSeq & HOMOZYGOUS] <- "C"
  groups[!is.na(HOMOZYGOUS) & N == 1 & !NoSeq & !HOMOZYGOUS] <- "E"
  groups[!is.na(HOMOZYGOUS) & N == 1 & !NoSeq & HOMOZYGOUS] <- "A"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & HOMOZYGOUS] <- "B"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & N_unique == 1 & !HOMOZYGOUS] <- "F"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & N_unique == 2 & !HOMOZYGOUS] <- "D"
  return(groups)
}

rename_methods <- function(method_names)
{
  method_names[method_names == "KMEANS"] <- "k-mer clust"
  method_names[method_names == "TRUE_HAP"] <- "True phase"
  method_names[method_names == "ITER_POA"] <- "Multi POA"
  method_names[method_names == "fuzzy"] <- "RapidFuzz"
  method_names[method_names == "Molephase"] <- "A priori"
  return(method_names)
}


# --- I/O variables --- #
variant_file <- "data/variant_regions.tsv"
results_dir <- "results"
output_dir <- "/home/tim/Git/lr_consensus_paper_code/results/plots/"
numbers_dir <- "/home/tim/Git/lr_consensus_paper_code/results/numbers/"
dir.create(output_dir)

# --- Load and preprocess data --- #

# load variant info, esp. zygosity and type
variant_regions <- read_tsv(variant_file, col_names = FALSE)
colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "ZYGOSITY", "REGION")
variant_regions <- variant_regions %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = toupper(ZYGOSITY)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM")

if (!file.exists(paste0(results_dir, "/results.rds"))) {
  # load evaluation results for all methods
  results <- tibble()
  for (technology in c("ont", "pacbio"))
  {
    acc_dirs <- list.dirs(paste0(results_dir, "/", technology), recursive = FALSE)
    for (acc_dir in acc_dirs)
    {
      acc <- tail(strsplit(acc_dir, split = "/")[[1]], n = 1)
      cov_dirs <- list.dirs(acc_dir, full.names = TRUE, recursive = FALSE)
      for (cov_dir in cov_dirs) {
        aln_dirs <- list.dirs(cov_dir, full.names = TRUE, recursive = FALSE)
        # extract coverage
        cov <- as.integer(tail(strsplit(cov_dir, split = "/")[[1]], n = 1))
        if (is.na(cov)) {
          next
        }
        for (aln_dir in aln_dirs) {
          aln_params <- tail(strsplit(aln_dir, split = "/")[[1]], n = 1)
          # extract aln_param string
          eval_files <- list.files(
            aln_dir,
            full.names = TRUE,
            pattern = ".*evaluation_.*tsv"
          )
          
          for (ef in eval_files) {
            aln_type <- "local"
            if (grepl(ef, pattern = "global")) {
              if (grepl(ef, pattern = "semiglobal"))
              {
                aln_type <- "semiglobal"
              } else {
                aln_type <- "global"
              }
            }
            spanning <- ifelse(grepl(ef, pattern = "spanning") || (aln_type == "global"), "spanning", "partial")
            if (grepl(ef, pattern = "merge_v1"))
            {
              stop("Unexpected file name: merge_v1")
            }
            
      	    if (grepl(ef, pattern = "merge_v2"))
            {
              stop("Unexpected file name: merge_v2")
            }
            
            info <- strsplit(ef, split = "/")[[1]]
            info <- gsub(info[length(info)], pattern = ".tsv", replacement = "")
            info <- gsub(info, pattern = "_semiglobal", replacement = "")
            info <- gsub(info, pattern = "_global", replacement = "")
            info <- gsub(info, pattern = "_spanning", replacement = "")
            
            
            method <- gsub(info, pattern = "evaluation_", replacement = "")
            cat("\nMethod: ", method, "\n\n")
            
            results <- results |> rbind(
              read_tsv(ef, col_names = TRUE) |>
                mutate(VARIANT = factor(VARIANT, levels = variant_regions$VARIANT)) |>
                complete(VARIANT) |>
                mutate(
                  COVERAGE = cov, 
                  ALN_PARAMS = aln_params,
                  TECHNOLOGY = technology,
                  ACCURACY = acc,
                  SPANNING = spanning,
                  ALN_TYPE = aln_type,
                  METHOD = method
                )
            )
          }
        }
      }
    }
  }

  # no longer needed, but created at some point during testing
  results <- results %>%
    filter(
      METHOD != "abPOA_conservative", 
      METHOD != "fuzzy_v2", 
      METHOD != "fuzzy_hclust", 
      METHOD != "ITER_POA_ONT",
      METHOD != "fuzzy_v3",
      METHOD != "ITER_POA_INTERMEDIATE"
    ) %>%
    filter(
      ! (METHOD == "ITER_POA" & TECHNOLOGY == "ont"),
      ! (METHOD == "ITER_POA_PERMISSIVE" & TECHNOLOGY == "pacbio")
    ) %>%
    mutate(
      METHOD = ifelse(TECHNOLOGY == "ont" & METHOD == "ITER_POA_PERMISSIVE", "ITER_POA", METHOD)
    )
  
  # extract variant size from name, calculate alignment match rate and join TYPE and ZYGOSITY
  results <- results %>% 
    mutate(
      VARSIZE = sapply(
        as.character(VARIANT),
        function(x) {
          return(as.numeric(strsplit(x, split = "_")[[1]][length(strsplit(x, split = "_")[[1]])]))
        }
      )
    ) %>%
    mutate(MATCH_RATE = MATCHED_BASES / CONSENSUS_LENGTH) %>%
    inner_join(variant_regions, by = c("VARIANT")) %>%
    mutate(COMPHET = grepl(TYPE, pattern = "COMPHET"))

  
  results <- results %>%
    mutate(TECHNOLOGY = toupper(TECHNOLOGY))
  results[which(results$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"

  results <- results %>%
    mutate(METHOD = rename_methods(METHOD))
  
  write_rds(results, paste0(results_dir, "/results.rds"))
} else {
  results <- read_rds(paste0(results_dir, "/results.rds"))
}

# Group results by variant zygosity, number and assignment of created consensus sequences
if (!file.exists(paste0(results_dir, "/grouped_variants.rds"))) {
  grouped_variants <- results %>%
    group_by(
      VARIANT, METHOD, ALN_PARAMS, COVERAGE, TECHNOLOGY, ACCURACY, VARSIZE, HOMOZYGOUS, ALN_TYPE, SPANNING
    ) %>%
    summarize(
      N = n(),
      N_unique = length(unique(ALLELE_MATCH)),
      NoSeq = n() == 1 & any(is.na(ALLELE_MATCH)),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(GROUP = determine_group(HOMOZYGOUS, N, N_unique, NoSeq))
  write_rds(grouped_variants, paste0(results_dir, "/grouped_variants.rds"))
} else {
  grouped_variants <- read_rds(paste0(results_dir, "/grouped_variants.rds"))
}

# Filter results down to true positive consensus sequences for Error Rate plots
# Only one sequence for a variant -> Keep that entry
# Two sequences for a homozygous variant: Keep the better error rate
# Two sequences for a heterozygous variant:
#   - Distinct: keep both entries
#   - Same allele match: keep the better entry
if (!file.exists(paste0(results_dir, "/true_positives.rds")))
{
  true_positives <- results %>%
    filter(!is.na(MATCH_RATE)) %>%
    group_by(VARIANT, METHOD, ALN_PARAMS, COVERAGE, TYPE, VARSIZE, HOMOZYGOUS, TECHNOLOGY, ACCURACY, ALN_TYPE, SPANNING, ALLELE_MATCH) |>
    summarize(
      ErrorRate = min(1 - MATCH_RATE), 
      .groups = "keep"
      ) |>
    group_by(VARIANT, METHOD, ALN_PARAMS, COVERAGE, TYPE, VARSIZE, HOMOZYGOUS, TECHNOLOGY, ACCURACY, ALN_TYPE, SPANNING) |>
    summarize(
      ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
      .groups = "keep"
    )
  write_rds(true_positives, paste0(results_dir, "/true_positives.rds"))
} else {
  true_positives <- read_rds(paste0(results_dir, "/true_positives.rds"))
}

# --- Plotting --- #

# --- Error Rate by Coverage --- #

true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "partial")|>
  group_by(METHOD, TECHNOLOGY) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |> write_tsv(paste0(numbers_dir, "/error_rates_tp_overall.tsv"), col_names = TRUE)

true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "spanning")|>
  group_by(METHOD, TECHNOLOGY) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |> write_tsv(paste0(numbers_dir, "/error_rates_tp_spanning.tsv"), col_names = TRUE)

true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "partial")|>
  group_by(METHOD, TECHNOLOGY, COVERAGE, VARSIZE) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |> 
  mutate(COVERAGE = 2*COVERAGE) |> 
  write_tsv(paste0(numbers_dir, "/error_rates_tp_size_cov.tsv"), col_names = TRUE)

# --- Error Rate averaged across all coverages and SV sizes --- #
p_er <- true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "partial", METHOD != "ITER_POA_PERMISSIVE") |>
  group_by(METHOD, TECHNOLOGY) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot(aes(x = TECHNOLOGY)) +
  geom_bar(
    aes(fill = factor(METHOD, levels = method_order), y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  geom_errorbar(
    aes(ymin = AvgER - seER, ymax = AvgER + seER, col = factor(METHOD, levels = method_order)), 
    position = "dodge"
  ) +
  guides(
    fill = guide_legend(title = "Method"),
    col = guide_none()
  ) +
  custom_theme + 
  ylab("Error Rate [%]") +
  scale_fill_manual(values = my_colors, limits = method_order) +
  scale_color_manual(values = my_colors, limits = method_order) +
  theme(
    axis.title.x = element_blank()
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed", alpha = 0.5) +
  theme(
    panel.grid.major.x = element_blank()
  )

pdf(paste0(output_dir, "/avg_er.pdf"), width = 10, height = 3)
plot(p_er)
dev.off()

# --- Number of variants per evaluation group --- #
p_counts <- grouped_variants |>
  filter(
    COVERAGE == 15, 
    ALN_PARAMS == "1_-1_-1_-1",
    ACCURACY == "0.99", 
    METHOD != "ITER_POA_PERMISSIVE",
    ALN_TYPE == "local", SPANNING == "partial"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  complete(GROUP, METHOD, TECHNOLOGY, fill = list(N = 0)) %>%
  ggplot(
    aes(
      x = as.factor(GROUP), 
      y = N, 
      fill = factor(METHOD, levels = method_order)
      )
  ) +
  facet_grid(rows = vars(TECHNOLOGY)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  scale_x_discrete(
    limits = factor(c("A", "B", "C", "D", "E", "F", "G")), 
    labels = group_names
  ) +
  ylab("# variants") +
  scale_fill_manual(limits = method_order, values = my_colors) +
  guides(fill = guide_legend(title = element_blank())) +
  geom_vline(
    xintercept = c(1.5, 2.5, 4.5, 5.5, 6.5), 
    linetype = "dashed", linewidth = 0.5, alpha = 0.6
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 1) +
  custom_theme +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

grouped_variants |>
  filter(
    COVERAGE == 15, 
    ALN_PARAMS == "1_-1_-1_-1",
    ACCURACY == "0.99", 
    METHOD != "ITER_POA_PERMISSIVE",
    ALN_TYPE == "local", SPANNING == "partial"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>% write_tsv(paste0(numbers_dir, "/group_counts.tsv"), col_names = TRUE)

pdf(paste0(output_dir, "/group_counts.pdf"), width = 25, height = 6)
plot(p_counts)
dev.off()


# --- Precision-Recall curves --- #
pr_values <- grouped_variants %>%
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99") %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(
    TP = Group_A + Group_B + 2*Group_D + Group_E + Group_F,
    FN = Group_C + Group_E + Group_F + 2*Group_G,
    FP = Group_B + Group_F,
    TN = 0,
    .groups = "keep"
  ) %>%
  mutate(
    Recall = TP / (TP + FN),
    Precision = TP / (TP + FP)
  ) %>%
  mutate(
    F1 = 2 * Precision * Recall / (Precision + Recall)
  ) %>% 
  mutate(COVERAGE = 2 * COVERAGE) %>%
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal"))

pr_values %>%
	write_tsv(paste0(numbers_dir, "/pr_values.tsv"), col_names = TRUE)

p_pr <- pr_values %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme +
  scale_y_continuous(breaks = c(0.9, 1.0)) #+
  # scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.9, 1.0))


p_pr_zoom_ls <- pr_values %>%
  filter(SPANNING == "spanning", TECHNOLOGY != "ONT", ALN_TYPE == "local") %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  geom_path(linewidth = 2) +
  geom_point(size = 5) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_none(),
    shape = guide_none()
  ) +
  custom_theme +
  theme(
    axis.title = element_blank()
  )

p_pr_zoom_gs <- pr_values %>%
  filter(SPANNING == "spanning", TECHNOLOGY != "ONT", ALN_TYPE == "global") %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  geom_path(linewidth = 2) +
  geom_point(size = 5) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_none(),
    shape = guide_none()
  ) +
  custom_theme +
  theme(
    axis.title = element_blank()
  )

pdf(paste0(output_dir, "/pr_plot_all.pdf"), width = 20, height = 6)
plot(p_pr)
dev.off()

pdf(paste0(output_dir, "/pr_plot_spanning_zoom.pdf"), width = 8, height = 5)
plot(p_pr_zoom_gs)
plot(p_pr_zoom_ls)
dev.off()


# --- Precision-Recall curves heterozygous--- #

pr_values_het <- grouped_variants %>%
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", !HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(
    TP = 2*Group_D + Group_E + Group_F,
    FN = Group_E + Group_F + 2*Group_G,
    FP = Group_F,
    TN = 0,
    .groups = "keep"
  ) %>%
  mutate(
    Recall = TP / (TP + FN),
    Precision = TP / (TP + FP)
  ) %>%
  mutate(
    F1 = 2 * Precision * Recall / (Precision + Recall)
  ) %>% 
  mutate(COVERAGE = 2 * COVERAGE) %>%
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal"))

pr_values_het %>%
  write_tsv(paste0(numbers_dir, "/pr_values_heteroygous.tsv"), col_names = TRUE)

p_pr_het <- pr_values_het %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme +
  ylim(0.9,1)

p_pr_zoom_ls_het <- pr_values_het %>%
  filter(SPANNING == "spanning", TECHNOLOGY != "ONT", ALN_TYPE == "local") %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  geom_path(linewidth = 2) +
  geom_point(size = 5) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_none(),
    shape = guide_none()
  ) +
  custom_theme +
  theme(
    axis.title = element_blank()
  )

p_pr_zoom_gs_het <- pr_values_het %>%
  filter(SPANNING == "spanning", TECHNOLOGY != "ONT", ALN_TYPE == "global") %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  geom_path(linewidth = 2) +
  geom_point(size = 5) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_none(),
    shape = guide_none()
  ) +
  custom_theme +
  theme(
    axis.title = element_blank()
  )

pdf(paste0(output_dir, "/pr_plot_all_het.pdf"), width = 20, height = 6)
plot(p_pr_het)
dev.off()

pdf(paste0(output_dir, "/pr_plot_spanning_zoom_het.pdf"), width = 8, height = 5)
plot(p_pr_zoom_gs_het)
plot(p_pr_zoom_ls_het)
dev.off()

# --- F1 score by variant type and size -- #



# --- Load time data and create a plot / write values --- #
stats <- tibble()
for (cov in c(10, 15, 30))
{
  filenames <- list.files(paste0("data/time_data/", cov), pattern = "*.tsv", recursive = FALSE)
  
  for (fn in filenames)
  {
    f <- paste0("data/time_data/", cov, "/", fn)
    
    fn <- gsub(fn, pattern = ".tsv", replacement = "")
    fn <- gsub(fn, pattern = "consensus_stats_", replacement = "")
    
    info <- strsplit(fn, split = "_")[[1]]
    
    aln_params <- paste0(info[(length(info) - 3) : length(info)], collapse = "_")
    accuracy <- info[(length(info) - 4)]
    technology <- info[(length(info) - 5)]
    method <- info[1:(length(info)-5)]
    
    temp_df <- read_tsv(
      f, col_names = c("VARIANT", "METHOD", "TIME", "FLAG", "SEQUENCE_COUNTS"), 
      col_types = list(VARIANT = "character", METHOD = "character", FLAG = "double", TIME = "double", SEQUENCE_COUNTS = "character")
    )
    stats <- stats |>
      rbind(
        temp_df |> mutate(
          ALN_PARAMS = aln_params,
          ACCURACY = accuracy,
          TECHNOLOGY = technology,
          COVERAGE = 2 * cov
        )
      )
  }
}

stats <- stats %>%
  mutate(
    VARSIZE = sapply(VARIANT, function(v){
      return(as.numeric(tail(strsplit(v, split = "_")[[1]], n = 1)))
    })
  )
stats <- stats %>%
  mutate(METHOD = rename_methods(METHOD))
# --- Write to disk --- #
stats %>%
  group_by(METHOD, ALN_PARAMS, TECHNOLOGY, ACCURACY, COVERAGE) %>%
  summarize(
    AvgTime = mean(TIME),
    seTime = sd(TIME) / sqrt(n()),
    .groups = "keep"
  ) %>% write_tsv(paste0(numbers_dir, "/time_averages_16000_ont_99_1_-1_-1_-1.tsv"))

p_time <- stats %>%
  group_by(METHOD, ALN_PARAMS, TECHNOLOGY, ACCURACY, COVERAGE) %>%
  summarize(
    AvgTime = mean(TIME),
    seTime = sd(TIME) / sqrt(n()),
    .groups = "keep"
  ) %>%
  ggplot(aes(x = COVERAGE, col = METHOD)) +
  geom_path(aes(y = AvgTime), linewidth = 2) +
  custom_theme +
  scale_color_manual(breaks = method_order, values = my_colors) +
  scale_x_continuous(breaks = c(20, 30, 60)) +
  xlab("Coverage") +
  ylab("Time [ms]")

pdf(paste0(output_dir, "/time_plot.pdf"), width = 12, height = 6)
plot(p_time)
dev.off()
