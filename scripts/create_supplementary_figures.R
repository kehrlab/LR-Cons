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
  if (any(method_names == "ITER_POA_PERMISSIVE"))
  {
    method_names[method_names == "Multi POA"] <- "Multi POA (0.96)"
    method_names[method_names == "ITER_POA_PERMISSIVE"] <- "Multi POA (0.9)"
  }
  return(method_names)
}


# --- I/O variables --- #
variant_file <- "data/variant_regions.tsv"
results_dir <- "results"
output_dir <- "results/supplementary_plots/"

dir.create(output_dir)


# load variant info, esp. zygosity and type
variant_regions <- read_tsv(variant_file, col_names = FALSE)
colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "ZYGOSITY", "REGION")
variant_regions <- variant_regions %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = toupper(ZYGOSITY)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM")


results <- read_rds(paste0(results_dir, "/results.rds"))
grouped_variants <- read_rds(paste0(results_dir, "/grouped_variants.rds"))
true_positives <- read_rds(paste0(results_dir, "/true_positives.rds"))

# --- Create alternative versions of paper plots --- #
pdf(paste0(output_dir, "/counts_coverages_local.pdf"), width = 25, height = 10)
for (p in c("1_-1_-1_-1", "1_-3_-1_-1", "1_-2_-3_-1"))
{
  p_counts_all <- grouped_variants |>
    filter(
      ALN_PARAMS == p,
      ACCURACY == "0.99",
      ALN_TYPE == "local", SPANNING == "partial"
    ) |>
    group_by(GROUP, METHOD, TECHNOLOGY, COVERAGE) %>%
    summarize(
      N = n(),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    complete(GROUP, METHOD, TECHNOLOGY, COVERAGE, fill = list(N = 0)) %>%
    ggplot(
      aes(
        x = as.factor(GROUP), 
        y = N, 
        fill = factor(METHOD, levels = method_order)
      )
    ) +
    facet_grid(rows = vars(TECHNOLOGY), cols = vars(2*COVERAGE)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    scale_x_discrete(
      limits = factor(c("A", "B", "C", "D", "E", "F", "G"))
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
      axis.title.x = element_blank()
    )
  plot(p_counts_all)
}
dev.off()


p_counts_global <- grouped_variants |>
  filter(
    COVERAGE == 15,
    ACCURACY == "0.99",
    ALN_TYPE == "global", SPANNING == "spanning"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, ALN_PARAMS) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  complete(GROUP, METHOD, TECHNOLOGY, ALN_PARAMS, fill = list(N = 0)) %>%
  ggplot(
    aes(
      x = as.factor(GROUP), 
      y = N, 
      fill = factor(METHOD, levels = method_order)
    )
  ) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(ALN_PARAMS)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  scale_x_discrete(
    limits = factor(c("A", "B", "C", "D", "E", "F", "G"))
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
    axis.title.x = element_blank()
  )

pdf(paste0(output_dir, "/counts_coverages_global_15.pdf"), width = 25, height = 10)
plot(p_counts_global)
dev.off()


p_counts_span <- grouped_variants |>
  filter(
    ALN_PARAMS == "1_-1_-1_-1",
    ACCURACY == "0.99",
    ALN_TYPE == "local", SPANNING == "spanning"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, COVERAGE) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  complete(GROUP, METHOD, TECHNOLOGY, COVERAGE, fill = list(N = 0)) %>%
  ggplot(
    aes(
      x = as.factor(GROUP), 
      y = N, 
      fill = factor(METHOD, levels = method_order)
    )
  ) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(2 * COVERAGE)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  scale_x_discrete(
    limits = factor(c("A", "B", "C", "D", "E", "F", "G"))
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
    axis.title.x = element_blank()
  )
pdf(paste0(output_dir, "/counts_coverages_spanning.pdf"), width = 25, height = 10)
plot(p_counts_span)
dev.off()


# --- Create an F1 score overview --- #
pr_values <- grouped_variants %>%
  filter(ACCURACY == "0.99") %>%
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

p_f1_overview <- pr_values %>%
  ggplot(aes(y = METHOD, x = F1)) +
  geom_bar(
    aes(
      fill = factor(interaction(ALN_TYPE, SPANNING, ALN_PARAMS, sep = " - "), levels = c(
        "local - partial - 1_-1_-1_-1",
        "local - partial - 1_-3_-1_-1",
        "local - partial - 1_-2_-3_-1",
        "semiglobal - partial - 1_-1_-1_-1",
        "local - spanning - 1_-1_-1_-1",
        "global - spanning - 1_-1_-1_-1",
        "global - spanning - 1_-3_-10_-1"
      ))
      ), 
    position = "dodge", stat = "identity") +
  facet_grid(cols = vars(TECHNOLOGY), rows = vars(COVERAGE)) +
  custom_theme + 
  scale_fill_manual(
    limits = c(
      "local - partial - 1_-1_-1_-1",
      "local - partial - 1_-3_-1_-1",
      "local - partial - 1_-2_-3_-1",
      "semiglobal - partial - 1_-1_-1_-1",
      "local - spanning - 1_-1_-1_-1",
      "global - spanning - 1_-1_-1_-1",
      "global - spanning - 1_-3_-10_-1"
    ),
    values = my_colors
  ) +
  guides(
    fill = guide_legend(title = "POA details")
  ) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  ) +
  scale_y_discrete(limits = method_order)

pdf(paste0(output_dir, "/f1_overview.pdf"), width = 20, height = 20)
plot(p_f1_overview)
dev.off()


# --- Create other versions of Fig. 3 --- #

# Homozygous only
pr_values_hom <- grouped_variants %>%
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", HOMOZYGOUS) %>% 
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS, ALN_TYPE, SPANNING, COVERAGE) %>%
  summarize(
    TP = Group_A + Group_B,
    FN = Group_C,
    FP = Group_B,
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


p_pr_hom <- pr_values_hom %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme

pdf(paste0(output_dir, "/pr_curves_hom_1_-1_-1_-1.pdf"), width = 16, height = 9)
plot(p_pr_hom)
dev.off()

# Hom & Het
pr_values_all <- grouped_variants %>%
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


p_pr_all <- pr_values_all %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme

pdf(paste0(output_dir, "/pr_curves_all.pdf"), width = 16, height = 9)
plot(p_pr_all)
dev.off()

# global
pr_values_global <- grouped_variants %>%
  filter(ACCURACY == "0.99", ALN_TYPE == "global") %>% 
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


p_pr_global <- pr_values_global %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, ALN_PARAMS, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme

pdf(paste0(output_dir, "/pr_curves_global.pdf"), width = 16, height = 9)
plot(p_pr_global)
dev.off()

# --- Fig 3. w
pr_values_aln_params <- grouped_variants %>%
  filter(ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "partial") %>% 
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_PARAMS, COVERAGE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS, COVERAGE) %>%
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
  mutate(COVERAGE = 2 * COVERAGE) 


p_pr_aln <- pr_values_aln_params %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_PARAMS, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme +
  xlim(0.3, 1)

pdf(paste0(output_dir, "/pr_aln_params.pdf"), width = 16, height = 9)
plot(p_pr_aln)
dev.off()


# --- create error rate plot per variant size, per coverage and overall for local and spanning --- # 
p_er_set <- true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local")|>
  group_by(METHOD, TECHNOLOGY, SPANNING) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |> 
  ggplot(aes(x = as.factor(0))) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(SPANNING)) +
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
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(
    panel.grid.major.x = element_blank()
  )


p_er_coverage <- true_positives |>
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99", ALN_TYPE == "local")|>
  group_by(METHOD, COVERAGE, TECHNOLOGY, SPANNING) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot(aes(x = METHOD)) +
  geom_bar(
    aes(fill = factor(2* COVERAGE, levels = c(20, 30, 60)), y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  geom_errorbar(
    aes(ymin = AvgER - seER, ymax = AvgER + seER, col = factor(2*COVERAGE, levels = c(20, 30, 60))), 
    position = "dodge"
  ) +
  guides(
    fill = guide_legend(title = "Coverage"),
    col = guide_none()
  ) +
  custom_theme + 
  rotate_text_x +
  ylab("Error Rate [%]") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = method_order) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(SPANNING)) +
  theme(
    panel.grid.major.x = element_blank()
  )

p_er_aln <- true_positives |>
  filter(ACCURACY == "0.99", ALN_TYPE == "local", SPANNING == "partial")|>
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot(aes(x = METHOD)) +
  geom_bar(
    aes(fill = ALN_PARAMS, y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  geom_errorbar(
    aes(ymin = AvgER - seER, ymax = AvgER + seER, col = ALN_PARAMS), 
    position = "dodge"
  ) +
  guides(
    fill = guide_legend(title = "POA parameters"),
    col = guide_none()
  ) +
  custom_theme + 
  rotate_text_x +
  ylab("Error Rate [%]") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = method_order) +
  facet_wrap(~TECHNOLOGY, ncol = 1) +
  theme(
    panel.grid.major.x = element_blank()
  )

p_er_global <- true_positives |>
  filter(ACCURACY == "0.99", ALN_TYPE == "global")|>
  group_by(METHOD, TECHNOLOGY, ALN_PARAMS) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot(aes(x = METHOD)) +
  geom_bar(
    aes(fill = ALN_PARAMS, y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  geom_errorbar(
    aes(ymin = AvgER - seER, ymax = AvgER + seER, col = ALN_PARAMS), 
    position = "dodge"
  ) +
  guides(
    fill = guide_legend(title = "Coverage"),
    col = guide_none()
  ) +
  custom_theme + 
  rotate_text_x +
  ylab("Error Rate [%]") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = method_order) +
  facet_wrap(~TECHNOLOGY, ncol = 1) +
  theme(
    panel.grid.major.x = element_blank()
  )
plot(p_er_global)

pdf(paste0(output_dir, "/error_rates.pdf"), width = 16, height = 9)
plot(p_er_set)
plot(p_er_coverage)
plot(p_er_aln)
plot(p_er_global)
dev.off()

# --- F1 score by variant type --- #
p1_t_f1 <- grouped_variants %>%
  filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == "0.99") %>%
  inner_join(variant_regions, by = "VARIANT") %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, TYPE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, TYPE) %>%
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
  filter((ALN_TYPE == "local" & SPANNING == "partial") | (ALN_TYPE == "global" & SPANNING == "spanning")) %>%
  ggplot(aes(
    x = METHOD, y = TYPE
  )) +
  geom_point(aes(size = F1, col = F1)) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  custom_theme +
  scale_x_discrete(limits = method_order) +
  rotate_text_x +
  scale_size_continuous(transform = "log", range = c(1, 10)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_gradient2(
    low = my_colors[7],
    mid = my_colors[4],
    high = my_colors[8],
    midpoint = 0.9
  ) +
  scale_y_discrete(
    limits = rev(c("DEL", "INV", "DUP", "COMPHET_DEL", "SHIFTED_COMPHET")), 
    labels = rev(c("Del", "Inv", "Dup", "Het_Del", "Het_Shift"))
    ) +
  guides(
    size = guide_none()
  )

pdf(paste0(output_dir, "/f1_by_type.pdf"), width = 12, height = 5)
plot(p1_t_f1)
dev.off()

# --- Sequence assignments --- #

# --- Function for count evaluation --- #
assigned_fraction <- function(SEQUENCE_COUNTS, TOTAL_SEQS) {
  all_counts <- strsplit(SEQUENCE_COUNTS, split = ",")
  main_counts <- sapply(all_counts, function(x) {
    return(sum(as.numeric(x[1:min(2, length(x)-1)])))
  })
  return(main_counts / TOTAL_SEQS)
}

unassigned_fraction <- function(SEQUENCE_COUNTS, TOTAL_SEQS) {
  all_counts <- strsplit(SEQUENCE_COUNTS, split = ",")
  all_counts <- sapply(all_counts, function(x) {
    return(sum(as.numeric(x[-length(x)])))
  })
  return(1 - (all_counts / TOTAL_SEQS))
}

assigned_distribution <- function(SEQUENCE_COUNTS) {
  all_counts <- strsplit(SEQUENCE_COUNTS, split = ",")
  y <- sapply(all_counts, function(x) {
    return(max(as.numeric(x[1:min(2, length(x)-1)])) / sum(as.numeric(x[1:min(2, length(x)-1)])))
  })
  return(y)
}

# --- ONT --- #
stats <- tibble()
for (cov in c(10, 15, 30))
{
  for (method in c("ITER_POA", "ITER_POA_PERMISSIVE", "fuzzy", "KMEANS", "abPOA", "Molephase", "TRUE_HAP"))
  {
    f <- paste0(
      "data/", cov, "/consensus_stats_merge_v2_", method, "_ont_0.99_1_-1_-1_-1.tsv"
    )
    temp_df <- read_tsv(
      f, 
      col_names = c("VARIANT", "METHOD", "TIME", "STATUS", "SEQUENCE_COUNTS"),
      col_types = list(VARIANT = "character", METHOD = "character", TIME = "double", STATUS = "double", SEQUENCE_COUNTS = "character")
    ) %>%
      mutate(COVERAGE = 2 * cov)
    stats <- stats |> rbind(temp_df)
  }
}
stats %>% filter(is.na(SEQUENCE_COUNTS)) %>% nrow
if (any(is.na(stats$SEQUENCE_COUNTS)))
{
  stats[which(is.na(stats$SEQUENCE_COUNTS)), "SEQUENCE_COUNTS"] <- "0,0,1"
}
stats <- stats |>
  mutate(
    TOTAL_SEQS = sapply(strsplit(SEQUENCE_COUNTS, split = ","), function(x) { return(as.numeric(x[length(x)])) })
  ) |>
  mutate(
    BALANCE = assigned_distribution(SEQUENCE_COUNTS),
    UNASSIGNED = unassigned_fraction(SEQUENCE_COUNTS, TOTAL_SEQS),
    ASSIGNED = assigned_fraction(SEQUENCE_COUNTS, TOTAL_SEQS),
    SECOND_N = sapply(strsplit(SEQUENCE_COUNTS, split = ","), function(x) { return(min(as.numeric(x[1:min(2, length(x) - 1)]))) })
  )
if (any(is.na(stats$BALANCE)))
{
  stats[which(is.na(stats$BALANCE)), "BALANCE"] <- 0.5
}


stats <- stats %>%
  inner_join(variant_regions, by = c("VARIANT")) %>%
  filter(!HOMOZYGOUS)

stats <- stats %>%
  mutate(METHOD = rename_methods(METHOD))

stats <- stats %>%
  mutate(METHOD = factor(METHOD, levels = c(
    "abPOA", "Multi POA (0.9)", "Multi POA (0.96)", "RapidFuzz", "k-mer clust", "A priori", "True phase"
  )))

p <- stats |>
  ggplot(aes(x = ASSIGNED, y = BALANCE)) +
  geom_point(aes(col = ifelse(SECOND_N >= 3, "2", "1")), alpha = 0.6) +
  geom_density2d(contour_var = "ndensity", col = "black") +
  custom_theme +
  scale_color_manual(values = my_colors[c(9,8)]) +
  facet_grid(cols = vars(METHOD), rows = vars(COVERAGE)) +
  theme(panel.spacing.x = unit(20, "pt")) +
  guides(col = guide_legend(title = "# Clusters"))

pdf(paste0(output_dir, "/read_assignments_ont.pdf"), width = 25, height =10)
plot(p)
dev.off()


# --- PacBio HiFi --- #
stats <- tibble()
for (cov in c(10, 15, 30))
{
  for (method in c("ITER_POA", "ITER_POA_PERMISSIVE", "fuzzy", "KMEANS", "abPOA", "Molephase", "TRUE_HAP"))
  {
    f <- paste0(
      "data/", cov, "/consensus_stats_merge_v2_", method, "_pacbio_0.99_1_-1_-1_-1.tsv"
    )
    temp_df <- read_tsv(
      f, 
      col_names = c("VARIANT", "METHOD", "TIME", "STATUS", "SEQUENCE_COUNTS"),
      col_types = list(VARIANT = "character", METHOD = "character", TIME = "double", STATUS = "double", SEQUENCE_COUNTS = "character")
    ) %>%
      mutate(COVERAGE = 2 * cov)
    stats <- stats |> rbind(temp_df)
  }
}
stats %>% filter(is.na(SEQUENCE_COUNTS)) %>% nrow
if (any(is.na(stats$SEQUENCE_COUNTS)))
{
  stats[which(is.na(stats$SEQUENCE_COUNTS)), "SEQUENCE_COUNTS"] <- "0,0,1"
}
stats <- stats |>
  mutate(
    TOTAL_SEQS = sapply(strsplit(SEQUENCE_COUNTS, split = ","), function(x) { return(as.numeric(x[length(x)])) })
  ) |>
  mutate(
    BALANCE = assigned_distribution(SEQUENCE_COUNTS),
    UNASSIGNED = unassigned_fraction(SEQUENCE_COUNTS, TOTAL_SEQS),
    ASSIGNED = assigned_fraction(SEQUENCE_COUNTS, TOTAL_SEQS),
    SECOND_N = sapply(strsplit(SEQUENCE_COUNTS, split = ","), function(x) { return(min(as.numeric(x[1:min(2, length(x) - 1)]))) })
  )
if (any(is.na(stats$BALANCE)))
{
  stats[which(is.na(stats$BALANCE)), "BALANCE"] <- 0.5
}


stats <- stats %>%
  inner_join(variant_regions, by = c("VARIANT")) %>%
  filter(!HOMOZYGOUS)

stats <- stats %>%
  mutate(METHOD = rename_methods(METHOD))

stats <- stats %>%
  mutate(METHOD = factor(METHOD, levels = c(
    "abPOA", "Multi POA (0.9)", "Multi POA (0.96)", "RapidFuzz", "k-mer clust", "A priori", "True phase"
  )))

p <- stats |>
  ggplot(aes(x = ASSIGNED, y = BALANCE)) +
  geom_point(aes(col = ifelse(SECOND_N >= 3, "2", "1")), alpha = 0.6) +
  geom_density2d(contour_var = "ndensity", col = "black") +
  custom_theme +
  scale_color_manual(values = my_colors[c(9,8)]) +
  facet_grid(cols = vars(METHOD), rows = vars(COVERAGE)) +
  theme(panel.spacing.x = unit(20, "pt")) +
  guides(col = guide_legend(title = "# Clusters"))

pdf(paste0(output_dir, "/read_assignments_pacbio.pdf"), width = 25, height =10)
plot(p)
dev.off()
