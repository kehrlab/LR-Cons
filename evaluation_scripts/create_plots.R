source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

### ---------------------------- Plotting ---------------------------------- ###

if (!dir.exists(output_dir))
{
  dir.create(output_dir)
}
if (!dir.exists(numbers_dir))
{
  dir.create(numbers_dir)
}

# Load variant regions
variant_regions <- read_tsv(variant_file_real, col_names = FALSE)
colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "GT", "REGION")
variant_regions <- variant_regions %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = zygosity_from_gt(GT)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM") %>% 
  mutate(SIZE = sapply(strsplit(VARIANT, split = "_"), function(x){return(as.integer(x[length(x) - 1]))})) %>%
  mutate(SIZE_BIN = (SIZE %/% 500) * 500)

# Load repeat annotations
variant_annot <- read_tsv(annotation_real)
variant_regions <- variant_regions %>%
  inner_join(variant_annot, by = "VARIANT")

# Load SV signal stats
region_stats <- read_tsv(stats_ont_real) %>%
  mutate(TECHNOLOGY = "ONT") %>%
  bind_rows(
    read_tsv(stats_pacbio_real) %>% 
      mutate(TECHNOLOGY = "PacBio HiFi")
  ) %>%
  mutate(SD_BIN = ifelse(SD_POS > 500, 500, (SD_POS %/% 50) * 50)) %>%
  select(VARIANT, TECHNOLOGY, SD_BIN, SD_POS)

# reduce to regions with sufficient coverage at 32.5x
no_low_cov <- read_tsv(variant_coverage_whitelist_old)

# variant_regions <- variant_regions %>% 
#   inner_join(no_low_cov, by = "VARIANT") %>% inner_join(no_overlap_500, by = "VARIANT")

# reduce to regions that do not overlap any other regions
no_overlap_500 <- read_tsv(variant_overlap_whitelist_500, col_names = c("VARIANT"))
variant_regions_500 <- variant_regions %>%
  inner_join(no_overlap_500, by = "VARIANT") %>%
  inner_join(no_low_cov, by = "VARIANT")
n_hom <- sum(variant_regions_500$HOMOZYGOUS)
n_het <- sum(!variant_regions_500$HOMOZYGOUS)
table(variant_regions_500$HOMOZYGOUS)


# group plot for Figure 2
grouped_variants_500 <- read_rds(paste0(results_dir_real, "/grouped_variants_real_500_v2.rds"))

p_counts <- grouped_variants_500 %>%
  mutate(TECHNOLOGY = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))) %>%
  filter(
    COVERAGE == 32.51,
    ALN_TYPE == "local", SPANNING == "partial",
    METHOD != "Multi POA v2"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, ORIGIN) %>%
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
      fill = factor(METHOD, levels = method_order_real)
    )
  ) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(ORIGIN)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  scale_x_discrete(
    limits = factor(c("A", "B", "C", "D", "E", "F", "G")), 
    labels = group_names
  ) +
  ylab("# variants") +
  scale_fill_manual(limits = method_order_real, values = my_colors) +
  guides(fill = guide_legend(title = element_blank())) +
  geom_vline(
    xintercept = c(1.5, 2.5, 4.5, 5.5, 6.5), 
    linetype = "dashed", linewidth = 0.5, alpha = 0.6
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 1) +
  custom_theme_poster +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

p_counts_ideal <- grouped_variants_500 %>%
  mutate(TECHNOLOGY = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))) %>%
  filter(ORIGIN == "real - 2025") %>%
  filter(
    COVERAGE == 32.51,
    ALN_TYPE == "global", SPANNING == "complete",
    METHOD != "Multi POA v2"
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
      fill = factor(METHOD, levels = method_order_real)
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
  scale_fill_manual(limits = method_order_real, values = my_colors) +
  guides(fill = guide_legend(title = element_blank())) +
  geom_vline(
    xintercept = c(1.5, 2.5, 4.5, 5.5, 6.5), 
    linetype = "dashed", linewidth = 0.5, alpha = 0.6
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 1) +
  custom_theme_poster +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

p_counts_comparison <- grouped_variants_500 %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA v2", "Multi POA (0.90)", METHOD)) %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA", "Multi POA (0.96)", METHOD)) %>%
  mutate(TECHNOLOGY = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))) %>%
  filter(
    COVERAGE == 32.51,
    ALN_TYPE == "global", SPANNING == "complete"#,
    # METHOD != "Multi POA v2"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, ORIGIN) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  complete(GROUP, METHOD, TECHNOLOGY, fill = list(N = 0)) %>%
  ggplot(
    aes(
      x = factor(METHOD), 
      y = N, 
      fill = ORIGIN
    )
  ) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(factor(GROUP, levels = c("A", "B", "C", "D", "E", "F", "G")))) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  ylab("# variants") + 
  custom_theme +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
  ) +
  scale_fill_manual(values = rev(my_colors)[c(1,4)])



pdf(paste0(output_dir, "/group_counts.pdf"), width = 25, height = 7)
plot(p_counts)
plot(p_counts_ideal) # this one goes in the paper
plot(p_counts_comparison)
dev.off()

# write corresponding numbers
grouped_variants_500 %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA v2", "Multi POA (0.90)", METHOD)) %>%
  mutate(TECHNOLOGY = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))) %>%
  filter(
    COVERAGE == 32.51,
    ALN_TYPE == "local", SPANNING == "partial"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, ORIGIN) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>% write_tsv(paste0(numbers_dir, "/group_counts_partial_local.tsv"))

grouped_variants_500 %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA v2", "Multi POA (0.90)", METHOD)) %>%
  mutate(TECHNOLOGY = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))) %>%
  filter(
    COVERAGE == 32.51,
    ALN_TYPE == "global", SPANNING == "complete",
    METHOD != "Multi POA v2"
  ) |>
  group_by(GROUP, METHOD, TECHNOLOGY, HOMOZYGOUS, ORIGIN) %>%
  summarize(
    N = n(),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  complete(GROUP, METHOD, TECHNOLOGY, fill = list(N = 0)) %>%
  mutate(N_total = ifelse(HOMOZYGOUS, n_hom, n_het)) %>%
  mutate(Fraction = (N / N_total) * 100) %>%
  write_tsv(paste0(numbers_dir, "/group_counts_spanning_global.tsv"))

# --- Create PR plot for real and simulated HiFi and ONT data ---
pr_values_real <- read_rds(paste0(results_dir_real, "/pr_values_real_v2.rds"))
pr_values_giab <- read_rds(paste0(results_dir_giab, "/pr_values_giab_full.rds"))
pr_values <- pr_values_real %>%
  bind_rows(pr_values_giab %>% mutate(ORIGIN = "pbsim"))

pr_values$SPANNING[pr_values$SPANNING == "complete"] <- "spanning"
pr_values$TECHNOLOGY[pr_values$TECHNOLOGY == "PacBio HiFi"] <- "HiFi"
pr_values <- pr_values %>%
  mutate(SPANNING = factor(SPANNING, levels = c("partial", "spanning"))) %>%
  mutate(ALN_TYPE = factor(ALN_TYPE, levels =c ("local", "semiglobal", "global")))

pr_values <- pr_values %>%
  filter(COVERAGE < 90)

pr_plot_full <- pr_values %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "), 
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim","ONT - 2025",  "ONT - 2020", "ONT - pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pr_plot_hifi_500 <- pr_values %>%
  filter(TECHNOLOGY == "HiFi") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "),
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  #xlim(0.45, 1) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pr_plot_hifi_1000 <- pr_values %>%
  filter(TECHNOLOGY == "HiFi") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "),
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim"))
  ) %>%
  filter(MARGIN == 1000) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pr_plot_hifi_500_zoomed <- pr_values %>%
  filter(TECHNOLOGY == "HiFi") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "),
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  #xlim(0.45, 1) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlim(75, 100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pdf(paste0(output_dir, "/pr_plot_hifi.pdf"), width = 24, height = 12)
plot(pr_plot_hifi_500)
plot(pr_plot_hifi_1000)
plot(pr_plot_hifi_500_zoomed)
dev.off()


pr_plot_ont_500 <- pr_values %>%
  filter(TECHNOLOGY == "ONT") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "),
      levels = c("ONT - 2025", "ONT - 2020", "ONT - pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlim(0,100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pr_plot_ont_1000 <- pr_values %>%
  filter(TECHNOLOGY == "ONT") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "), 
      levels = c("ONT - 2025", "ONT - 2020", "ONT - pbsim"))
  ) %>%
  filter(MARGIN == 1000) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlim(0,100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pr_plot_ont_500_zoomed <- pr_values %>%
  filter(TECHNOLOGY == "ONT") %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "),
      levels = c("ONT - 2025", "ONT - 2020", "ONT - pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlim(75, 100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pdf(paste0(output_dir, "/pr_plot_ont.pdf"), width = 24, height = 12)
plot(pr_plot_ont_500)
plot(pr_plot_ont_1000)
plot(pr_plot_ont_500_zoomed)
dev.off()

pdf(paste0(output_dir, "/pr_plot_full.pdf"), width = 24, height = 15)
plot(pr_plot_full)
dev.off()

pr_values$COVERAGE[pr_values$COVERAGE == 32.51] <- 30


pr_plot_paper <- pr_values %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " "), 
      levels = c("HiFi 2025", "HiFi 2020", "HiFi pbsim","ONT 2025",  "ONT 2020", "ONT pbsim"))
  ) %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(TECHNOLOGY == "ONT") %>%
  #mutate(ORIGIN =  ifelse(ORIGIN == "giab_sim", "sim", "real")) %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  #xlim(0.45, 1) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlim(0, 100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pdf(paste0(output_dir, "/pr_plot.pdf"), width = 25, height = 8)
plot(pr_plot_paper)
dev.off()

pr_plot_multi <- pr_values %>%
  filter(grepl(x = METHOD, pattern = "Multi POA")) %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA", "Multi POA (0.96)", "Multi POA (0.90)")) %>%
  filter((ALN_TYPE == "local" & SPANNING == "partial") | (ALN_TYPE == "global" & SPANNING == "spanning")) %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "), 
      # levels = c("HiFi - real - 2025", "ONT - real - 2025", "HiFi - real - 2020", "ONT - real - 2020", "HiFi - pbsim", "ONT - pbsim"))
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim","ONT - 2025",  "ONT - 2020", "ONT - pbsim"))
  ) %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, MARGIN, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = c("Multi POA (0.96)", "Multi POA (0.90)"), values = my_colors[c(3,9)]) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines"),
    axis.text.x = element_text(angle = 35, vjust = 0.9, hjust = 0.9)
  ) +
  xlim(0, 100) +
  xlab("Recall [%]") +
  ylab("Precision [%]")

pdf(paste0(output_dir, "/pr_plot_multipoa.pdf"), width = 24, height = 15)
plot(pr_plot_multi)
dev.off()

pr_values %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  write_tsv(paste0(numbers_dir, "/pr_values_plot.tsv")) 


pr_values %>%
  filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(ALN_TYPE == "global" | (ALN_TYPE == "local" & SPANNING == "partial")) %>%
  pivot_wider(id_cols = c("METHOD", "TECHNOLOGY", "ALN_TYPE", "SPANNING", "MARGIN", "ORIGIN"), names_from = COVERAGE, values_from = Recall) %>%
  mutate(Diff = (`30` - `15`) * 100) %>%
  write_tsv(paste0(numbers_dir, "/pr_value_coverage_increase_plot.tsv"))


pr_values %>%
  mutate(METHOD = ifelse(METHOD == "Multi POA v2", "Multi POA (0.90)", METHOD)) %>%
  write_tsv(paste0(numbers_dir, "/pr_values_all.tsv")) 


# --- plot margin effect on real data ---
margin_effect <- pr_values %>%
  filter(COVERAGE == 30) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(grepl("real", ORIGIN), SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_wider(id_cols = c("SPANNING", "ALN_TYPE", "TECHNOLOGY", "COVERAGE", "METHOD", "ORIGIN"), names_from = "MARGIN", names_prefix = "F1_", values_from = "F1") %>%
  mutate(F1_DIFF = F1_1000 - F1_500) %>%
  group_by(METHOD) %>%
  summarize(M = mean(F1_DIFF), S = sd(F1_DIFF), .groups = "keep")
method_order_margin <- margin_effect %>% arrange(M) %>% pull(METHOD)


p_margin <- pr_values %>%
  filter(COVERAGE == 30) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(! (METHOD == "Multi POA" & ORIGIN == "real - 2020" & TECHNOLOGY == "ONT")) %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  filter(grepl("real", ORIGIN), SPANNING == "spanning" & ALN_TYPE == "global") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot() +
  facet_grid(
    cols = vars(factor(Metric, levels = c("Recall", "Precision"))), 
    rows = vars(factor(interaction(TECHNOLOGY, YEAR, sep = " "), levels = c("HiFi 2020", "HiFi 2025", "ONT 2020", "ONT 2025"))), 
    scale = "free_x") +
  geom_point(aes(x = Value * 100, y = factor(METHOD, levels = method_order_margin), col = factor(MARGIN)), size = 4, alpha = 0.9) +
  scale_color_manual(values = rev(my_colors)[c(2, 4)]) +
  custom_theme_poster +
  guides(col = guide_legend(title = "Flanking region (bp)", position = "top", nrow = 1)) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("%") +
  xlim(80, 100)

p_margin_v2 <- pr_values %>%
  filter(COVERAGE == 30) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(! (METHOD == "Multi POA" & ORIGIN == "real - 2020" & TECHNOLOGY == "ONT")) %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  filter(grepl("real", ORIGIN), SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot() +
  facet_grid(
    cols = vars(factor(Metric, levels = c("Recall", "Precision"))), 
    rows = vars(factor(interaction(TECHNOLOGY, YEAR, sep = " "), levels = c("HiFi 2020", "HiFi 2025", "ONT 2020", "ONT 2025"))), 
    scale = "free_x") +
  geom_point(aes(x = Value*100, y = factor(METHOD, levels = method_order_margin), col = factor(MARGIN)), size = 4, alpha = 0.9) +
  scale_color_manual(values = rev(my_colors)[c(2, 4)]) +
  custom_theme_poster +
  guides(col = guide_legend(title = "Flanking region (bp)", position = "top", nrow = 1)) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("%") +
  xlim(70, 100)


pdf(paste0(output_dir, "/flanking_plot.pdf"), width = 16, height = 11)
plot(p_margin)
plot(p_margin_v2)
dev.off()

pr_values %>%
  filter(COVERAGE == 30) %>%
  filter(grepl("real", ORIGIN), SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  write_tsv(paste0(numbers_dir, "/pr_values_margin_effect_real.tsv"))

pr_values %>%
  filter(COVERAGE == 30) %>%
  filter(grepl("real", ORIGIN), SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>% 
  pivot_wider(id_cols = c("TECHNOLOGY", "Metric", "METHOD", "SPANNING", "ORIGIN"), names_from = "MARGIN", values_from = "Value", names_prefix = "Value_") %>%
  mutate(Difference = Value_1000 - Value_500) %>%
  write_tsv(paste0(numbers_dir, "/pr_value_difference_margin_real.tsv"))


### Create partial vs complete plot on all data

# load purely simulated data and compare

pr_values_all <- read_rds(paste0(results_dir_real, "/pr_values_real.rds")) %>%
  bind_rows(read_rds(paste0(results_dir_giab, "/pr_values_giab_full.rds")) %>% mutate(ORIGIN = "simulated - GIAB")) %>%
  bind_rows(
    read_rds(paste0(results_dir, "/pr_values.rds")) %>% mutate(MARGIN = 500) %>%
      filter(ALN_PARAMS == "1_-1_-1_-1", ACCURACY == 0.99) %>%
      mutate(ORIGIN = "simulated - uniform")
    )

pr_values_all$SPANNING[pr_values_all$SPANNING == "complete"] <- "spanning"
pr_values_all$ORIGIN <- factor(pr_values_all$ORIGIN, levels = c("real - 2020", "real - 2025", "simulated - GIAB", "simulated - uniform"))


partial_effect <- pr_values_all %>%
  filter(COVERAGE == 30 | COVERAGE == 32.51) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  filter(MARGIN == 500) %>%
  pivot_wider(id_cols = c("TECHNOLOGY", "COVERAGE", "METHOD", "ORIGIN"), names_from = "SPANNING", names_prefix = "F1_", values_from = "F1") %>%
  mutate(F1_DIFF = F1_partial - F1_spanning) %>%
  group_by(METHOD) %>%
  summarize(M = mean(F1_DIFF), S = sd(F1_DIFF), .groups = "keep")
method_order_partial <- partial_effect %>% arrange(M) %>% pull(METHOD)

p_partial <- pr_values_all %>%
  filter(COVERAGE == 30 | COVERAGE == 32.51) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  filter(MARGIN == 500) %>%
  # filter(ORIGIN != "simulated") %>%
  pivot_longer(cols = c("Recall", "Precision"), names_to = "Metric", values_to = "Value") %>%
  ggplot() +
  facet_grid(
    cols = vars(factor(Metric, levels = c("Recall", "Precision"))), 
    rows = vars(ORIGIN), 
    scale = "free_x") +
  geom_point(
    aes(
      x = Value * 100, y = factor(METHOD, levels = rev(c("True phase", "A priori", "RapidFuzz", "k-mer clust", "abPOA", "Multi POA"))), col = SPANNING, 
      shape = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))
      ), size = 4, alpha = 0.9) +
  custom_theme_poster +
  scale_color_manual(values = rev(my_colors)[c(1,6)]) +
  guides(
    col = guide_legend(title = "Input Set", position = "top", nrow = 1), 
    shape = guide_legend(title = "Technology", position = "top", nrow = 1)) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("Value [%]") 

p_partial_v2 <- pr_values_all %>%
  filter(COVERAGE == 30 | COVERAGE == 32.51) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(!(METHOD == "Multi POA" & TECHNOLOGY == "ONT" & ORIGIN == "real - 2020")) %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  filter(MARGIN == 500) %>%
  # filter(ORIGIN != "simulated") %>%
  pivot_longer(cols = c("Recall", "Precision"), names_to = "Metric", values_to = "Value") %>%
  ggplot() +
  facet_grid(
    cols = vars(factor(Metric, levels = c("Recall", "Precision"))), 
    rows = vars(ORIGIN), 
    scale = "free_x") +
  geom_point(
    aes(
      x = Value * 100, y = factor(METHOD, levels = rev(c("True phase", "A priori", "RapidFuzz", "k-mer clust", "abPOA", "Multi POA"))), col = SPANNING, 
      shape = factor(TECHNOLOGY, levels = c("PacBio HiFi", "ONT"))
    ), size = 4, alpha = 0.9) +
  custom_theme_poster +
  scale_color_manual(values = rev(my_colors)[c(1,6)]) +
  guides(
    col = guide_legend(title = "Input Set", position = "top", nrow = 1), 
    shape = guide_legend(title = "Technology", position = "top", nrow = 1)) +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("Value [%]") 


pdf(paste0(output_dir, "/partial_effect.pdf"), width = 16, height = 14)
plot(p_partial)
plot(p_partial_v2)
dev.off()

pr_values_all %>%
  filter(COVERAGE == 30 | COVERAGE == 32.51) %>%
  filter(METHOD != "Multi POA v2") %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  filter(MARGIN == 500) %>%
  pivot_longer(cols = c("Recall", "Precision"), names_to = "Metric", values_to = "Value") %>%
  write_tsv(paste0(numbers_dir, "/pr_values_partial_spanning.tsv"))


# --- create PR plot without filtering low-coverage regions

# pr_values_real_all <- read_rds(paste0(results_dir_real, "/pr_values_real_all_500.rds"))
# 
# pr_values_real_all$SPANNING[pr_values_real_all$SPANNING == "complete"] <- "spanning"
# pr_values_real_all$TECHNOLOGY[pr_values_real_all$TECHNOLOGY == "PacBio HiFi"] <- "HiFi"
# pr_values_real_all <- pr_values_real_all %>%
#   mutate(SPANNING = factor(SPANNING, levels = c("partial", "spanning"))) %>%
#   mutate(ALN_TYPE = factor(ALN_TYPE, levels =c ("local", "semiglobal", "global")))
# 
# pr_plot_all <- pr_values_real_all %>%
#   mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
#   filter(YEAR == "2025") %>%
#   mutate(
#     DATA = factor(
#       interaction(TECHNOLOGY, YEAR, sep = " - "), 
#       levels = c("HiFi - 2025", "HiFi - 2020", "ONT - 2025",  "ONT - 2020"))
#   ) %>%
#   filter(MARGIN == 500) %>%
#   filter(METHOD != "Multi POA v2") %>%
#   ggplot(aes(x = Recall, y = Precision, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
#   facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
#   geom_path(linewidth = 1) +
#   geom_point(size = 3) +
#   scale_color_manual(limits = method_order, values = my_colors) +
#   guides(
#     col = guide_legend(title = "Method", position = "top", nrow = 1),
#     shape = guide_legend(title = "Coverage", position = "top")
#   ) +
#   custom_theme_poster +
#   theme(
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
#     panel.spacing.y = unit(1, "lines")
#   ) +
#   xlim(0.75, 0.87)
# plot(pr_plot_all)
# 
# pr_values_real_all %>%
#   filter(MARGIN == 500) %>%
#   filter(METHOD != "Multi POA v2") %>%
#   write_tsv(paste0(numbers_dir, "/pr_values_with_low_cov.tsv"))

# --- create PR plot for homozygous variants
grouped_variants_giab <- read_rds(paste0(results_dir_giab, "/grouped_variants_giab_full.rds"))
pr_values_500_hom_real <- grouped_variants_500 %>%
  filter(HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
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
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal"))
pr_values_500_hom_giab <- grouped_variants_giab %>%
  filter(HOMOZYGOUS) %>%
  mutate(ORIGIN = "pbsim") %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
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
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal"))

pr_values_500_hom <- pr_values_500_hom_real %>%
  rbind(pr_values_500_hom_giab)

pr_plot_full_hom <- pr_values_500_hom %>%
  mutate(TECHNOLOGY = ifelse(TECHNOLOGY == "PacBio HiFi", "HiFi", "ONT")) %>%
  filter(COVERAGE < 90) %>%
  mutate(COVERAGE = ifelse(COVERAGE == 32.51, 30, COVERAGE)) %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  mutate(
    DATA = factor(
      interaction(TECHNOLOGY, YEAR, sep = " - "), 
      levels = c("HiFi - 2025", "HiFi - 2020", "HiFi - pbsim","ONT - 2025",  "ONT - 2020", "ONT - pbsim"))
  ) %>%
  # filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2") %>%
  ggplot(aes(x = Recall * 100, y = Precision * 100, group = METHOD, col = METHOD, shape = factor(COVERAGE))) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(DATA)) +
  geom_path(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(limits = method_order, values = my_colors) +
  guides(
    col = guide_legend(title = "Method", position = "top", nrow = 1),
    shape = guide_legend(title = "Coverage", position = "top")
  ) +
  custom_theme_poster +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6)),
    panel.spacing.y = unit(1, "lines")
  ) +
  xlab("Recall [%]") +
  ylab("Precision [%]")


pr_values_500_hom %>%
  mutate(TECHNOLOGY = ifelse(TECHNOLOGY == "PacBio HiFi", "HiFi", "ONT")) %>%
  filter(COVERAGE < 90) %>%
  mutate(COVERAGE = ifelse(COVERAGE == 32.51, 30, COVERAGE)) %>%
  mutate(YEAR = gsub(x = ORIGIN, pattern = "real - ", replacement = "")) %>%
  filter(YEAR == "2025") %>%
  # filter(MARGIN == 500) %>%
  filter(METHOD != "Multi POA v2", ALN_TYPE == "global", COVERAGE == 30) %>%
  write_tsv(paste0(numbers_dir, "/pr_values_for_2a_counts.tsv"))


pdf(paste0(output_dir, "/pr_plot_homozygous.pdf"), width = 24, height = 15)
plot(pr_plot_full_hom)
dev.off()

### create margin plot (for jitter 0, 50, 200)
pr_values_sim <- read_rds(paste0(results_dir_sim, "/pr_values_uncertainty.rds"))
pr_values_sim$SPANNING[pr_values_sim$SPANNING == "complete"] <- "spanning"
p_margin_sim_ont <- pr_values_sim %>%
  filter(COVERAGE == 30, MARGIN > 100, TECHNOLOGY == "ONT", METHOD != "Multi POA v2") %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  pivot_wider(id_cols = c("METHOD", "JITTER", "SPANNING", "ALN_TYPE", "COVERAGE", "TECHNOLOGY", "Metric"), names_from = "MARGIN", names_prefix = "Value", values_from = "Value") %>%
  mutate(Difference = Value1000 - Value500) %>%
  ggplot() +
  geom_bar(aes(y = METHOD, x = Difference * 100, fill = factor(JITTER)), stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors) +
  custom_theme +
  guides(fill = guide_legend(title = "Region uncertainty")) +
  facet_grid(rows = vars(Metric), cols = vars(SPANNING)) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  xlab("Difference [p.pt.]") +
  ggtitle("ONT")
p_margin_sim_hifi <- pr_values_sim %>%
  filter(COVERAGE == 30, MARGIN > 100, TECHNOLOGY != "ONT", METHOD != "Multi POA v2") %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  pivot_wider(id_cols = c("METHOD", "JITTER", "SPANNING", "ALN_TYPE", "COVERAGE", "TECHNOLOGY", "Metric"), names_from = "MARGIN", names_prefix = "Value", values_from = "Value") %>%
  mutate(Difference = Value1000 - Value500) %>%
  ggplot() +
  geom_bar(aes(y = METHOD, x = Difference * 100, fill = factor(JITTER)), stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors) +
  custom_theme +
  guides(fill = guide_legend(title = "Region uncertainty")) +
  facet_grid(rows = vars(Metric), cols = vars(SPANNING)) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  xlab("Difference [p.pt.]") +
  ggtitle("PacBio HiFi")

p_margin_recall <- pr_values_sim %>%
  filter(COVERAGE == 30, METHOD != "Multi POA v2", MARGIN > 100) %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global") %>%
  # pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  # pivot_wider(id_cols = c("METHOD", "JITTER", "SPANNING", "ALN_TYPE", "COVERAGE", "TECHNOLOGY", "Metric"), names_from = "MARGIN", names_prefix = "Value", values_from = "Value") %>%
  # mutate(Difference = Value1000 - Value500) %>%
  # filter(Metric == "Recall") %>%
  ggplot() +
  geom_point(aes(x = Recall*100, col = factor(MARGIN), y = METHOD), size = 3) +
  facet_grid(rows = vars(JITTER), cols = vars(TECHNOLOGY)) +
  custom_theme_poster +
  theme(
    axis.title.y = element_blank()
  ) +
  xlim(85, 100) +
  scale_y_discrete(limits = rev(c("True phase", "A priori", "RapidFuzz", "k-mer clust", "abPOA", "Multi POA"))) +
  scale_color_manual(breaks = c(500, 1000), values = rev(my_colors)[c(2,4)]) +
  guides(col = guide_legend(title = "Flanking region (bp)", position = "top")) +
  xlab("Recall [%]")

pdf(paste0(output_dir, "/recall_margin_sim_uncertainty.pdf"), width = 16, height = 11)
plot(p_margin_recall)
dev.off()

pdf(paste0(output_dir, "/pr_margin_sim_uncertainty.pdf"), width = 16, height = 8)
plot(p_margin_sim_ont)
plot(p_margin_sim_hifi)
dev.off()

pr_values_sim %>%
  filter(COVERAGE == 30, MARGIN > 100) %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING == "partial" & ALN_TYPE == "local") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  pivot_wider(id_cols = c("METHOD", "JITTER", "SPANNING", "ALN_TYPE", "COVERAGE", "TECHNOLOGY", "Metric"), names_from = "MARGIN", names_prefix = "Value", values_from = "Value") %>%
  mutate(Difference = Value1000 - Value500) %>%
  write_tsv(paste0(numbers_dir, "/pr_diff_margin_by_uncertainty.tsv"))


# --- create time plot

stats <- tibble()

filenames <- list.files(paste0("data/time_benchmark/complete/global/"), pattern = "consensus_stats_.*.tsv", recursive = FALSE)
for (fn in filenames)
{
  f <- paste0("data/time_benchmark/complete/global/", fn)

  fn <- gsub(fn, pattern = ".tsv", replacement = "")
  fn <- gsub(fn, pattern = "consensus_stats_", replacement = "")
  
  info <- strsplit(fn, split = "_")[[1]]
  
  margin <- as.numeric(info[length(info)])
  cov <- as.numeric(info[length(info) - 1])
  technology <- info[(length(info) - 2)]
  method <- paste(info[1:(length(info)-3)], collapse = "_")
  
  temp_df <- read_tsv(
    f, col_names = c("VARIANT", "METHOD", "TIME", "FLAG", "SEQUENCE_COUNTS"), 
    col_types = list(VARIANT = "character", METHOD = "character", FLAG = "double", TIME = "double", SEQUENCE_COUNTS = "character")
  )
  stats <- stats |>
    rbind(
      temp_df |> mutate(
        TECHNOLOGY = technology,
        COVERAGE = cov,
        MARGIN = margin, 
        ALN_TYPE = "global",
        SPANNING = "spanning"
      )
    )
}

filenames <- list.files(paste0("data/time_benchmark/partial/local/"), pattern = "consensus_stats_.*.tsv", recursive = FALSE)
for (fn in filenames)
{
  f <- paste0("data/time_benchmark/partial/local/", fn)
  
  fn <- gsub(fn, pattern = ".tsv", replacement = "")
  fn <- gsub(fn, pattern = "consensus_stats_", replacement = "")
  
  info <- strsplit(fn, split = "_")[[1]]
  
  margin <- as.numeric(info[length(info)])
  cov <- as.numeric(info[length(info) - 1])
  technology <- info[(length(info) - 2)]
  method <- paste(info[1:(length(info)-3)], collapse = "_")
  
  temp_df <- read_tsv(
    f, col_names = c("VARIANT", "METHOD", "TIME", "FLAG", "SEQUENCE_COUNTS"), 
    col_types = list(VARIANT = "character", METHOD = "character", FLAG = "double", TIME = "double", SEQUENCE_COUNTS = "character")
  )
  stats <- stats |>
    rbind(
      temp_df |> mutate(
        TECHNOLOGY = technology,
        COVERAGE = cov,
        MARGIN = margin,
        ALN_TYPE = "local",
        SPANNING = "partial"
      )
    )
}


# stats <- stats %>%
#   mutate(
#     VARSIZE = sapply(VARIANT, function(v){
#       return(as.numeric(tail(strsplit(v, split = "_")[[1]], n = 1)))
#     })
#   )
stats <- stats %>%
  mutate(METHOD = rename_methods(METHOD))
stats <- stats %>%
  filter(METHOD != "Multi POA v2") %>%
  mutate(TECHNOLOGY = ifelse(TECHNOLOGY == "ont", "ONT", "PacBio HiFi"))

# --- Write to disk --- #
stats %>%
  group_by(METHOD, TECHNOLOGY, MARGIN, COVERAGE, ALN_TYPE, SPANNING) %>%
  summarize(
    AvgTime = mean(TIME),
    seTime = sd(TIME) / sqrt(n()),
    .groups = "keep"
  ) %>% write_tsv(paste0(numbers_dir, "/time_averages.tsv"))

p_time <- stats %>%
  group_by(METHOD, TECHNOLOGY, MARGIN, COVERAGE, ALN_TYPE, SPANNING) %>%
  summarize(
    AvgTime = mean(TIME),
    seTime = sd(TIME) / sqrt(n()),
    .groups = "keep"
  ) %>%
  ggplot(aes(x = COVERAGE, col = METHOD)) +
  geom_path(aes(y = AvgTime), linewidth = 2) +
  facet_grid(rows = vars(MARGIN), cols = vars(interaction(TECHNOLOGY, SPANNING, sep = " - "))) +
  scale_color_manual(breaks = method_order, values = my_colors) +
  scale_x_continuous(breaks = c(20, 30, 60)) +
  xlab("Coverage") +
  ylab("Time [ms]") +
  custom_theme

pdf(paste0(output_dir, "/time_plot.pdf"), width = 16, height = 9)
plot(p_time)
dev.off()



# --- Missed variants


############# Supplementary plots ###################################

output_dir <- "results/plots_final/supplement"
numbers_dir <- "results/numbers_final/supplement"  
if (!dir.exists(output_dir))
{
  dir.create(output_dir)
}
if (!dir.exists(numbers_dir))
{
  dir.create(numbers_dir)
}

pr_values_kmeans <- read_rds("results/kmeans_simulation/pr_values_kmeans.rds")
pr_values_kmeans$SPANNING[pr_values_kmeans$SPANNING == "complete"] <- "spanning" 
p_pr <- pr_values_kmeans %>%
  filter(ACCURACY == 0.99) %>%
  ggplot(aes(x = Recall, y = Precision, group = METHOD)) +
  facet_grid(cols = vars(interaction(ALN_TYPE, SPANNING, sep = " - ")), rows = vars(TECHNOLOGY)) +
  geom_path(linewidth = 1) +
  geom_point(aes(col = factor(KMER_SIZE)),size = 3) +
  scale_color_manual(limits = factor(c(5, 7, 11, 21, 31)), values = my_colors) +
  guides(
    col = guide_legend(title = "K", position = "top", nrow = 1)
  ) +
  custom_theme +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = scales::alpha("grey92", .6))
  )
pdf(paste0(output_dir, "/pr_kmeans.pdf"), width = 16, height = 7)
plot(p_pr)
dev.off()


true_pos_kmeans <- read_rds("results/kmeans_simulation/true_positives.rds")
true_pos_kmeans$SPANNING[true_pos_kmeans$SPANNING == "complete"] <- "spanning"
p_er <- true_pos_kmeans |>
  group_by(TECHNOLOGY, SPANNING, KMER_SIZE) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot(aes(x = TECHNOLOGY)) +
  facet_grid(rows = vars(SPANNING)) +
  geom_bar(
    aes(fill = factor(KMER_SIZE, levels = c(5, 7, 11, 21, 31)), y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  geom_errorbar(
    aes(ymin = AvgER - seER, ymax = AvgER + seER, col = factor(KMER_SIZE, levels = c(5, 7, 11, 21, 31))), 
    position = "dodge"
  ) +
  guides(
    fill = guide_legend(title = "Method"),
    col = guide_none()
  ) +
  custom_theme + 
  ylab("Error Rate [%]") +
  scale_fill_manual(values = my_colors, limits = factor(c(5, 7, 11, 21, 31))) +
  scale_color_manual(values = my_colors, limits = factor(c(5, 7, 11, 21, 31))) +
  theme(
    axis.title.x = element_blank()
  ) +
  geom_vline(xintercept = 1.5, linetype = "dashed", alpha = 0.5) +
  theme(
    panel.grid.major.x = element_blank()
  )

pdf(paste0(output_dir, "/error_rate_kmeans.pdf"), width = 16, height = 7)
plot(p_er)
dev.off()

# real - spanning, local
true_positives_real <- read_rds(paste0(results_dir_real, "/true_positives_real_v2.rds"))

true_positives_real$SPANNING[true_positives_real$SPANNING == "complete"] <- "spanning"
true_positives_real <- true_positives_real %>%
  filter(COVERAGE == 32.51) %>%
  filter(SPANNING == "spanning" & ALN_TYPE == "global" | SPANNING  == "partial" & ALN_TYPE == "local") 
true_positives_real <- true_positives_real %>%
  filter(METHOD != "Multi POA v2")


p_er_real <- true_positives_real %>%
  ggplot(aes(x = METHOD, y = ErrorRate * 100, col = METHOD)) +
  facet_grid(cols = vars(interaction(TECHNOLOGY, ORIGIN, sep = " - ")), rows = vars(SPANNING)) +
  geom_violin() +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
  ) +
  ylab("Error Rate [%]")

p_er_real_v2 <- true_positives_real %>%
  group_by(METHOD, TECHNOLOGY, ORIGIN, SPANNING) %>%
  summarize(AvgER = median(ErrorRate * 100), .groups = "keep") %>%
  ggplot(aes(x = METHOD, y = AvgER, fill = METHOD)) +
  facet_grid(cols = vars(interaction(TECHNOLOGY, ORIGIN, sep = " - ")), rows = vars(SPANNING)) +
  geom_bar(stat = "identity") +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
  ) +
  ylab("Error Rate [%]") 

pdf(paste0(output_dir, "/error_rate_real.pdf"), width = 16, height= 7)
plot(p_er_real)
plot(p_er_real_v2)
dev.off()

error_rate_real <- true_positives_real %>%
  group_by(METHOD, SPANNING, TECHNOLOGY, ORIGIN) %>%
  summarize(AvgER = median(ErrorRate), .groups = "keep")
error_rate_real %>% write_tsv(paste0(numbers_dir, "/median_error_rates_real.tsv"))


## repeat context

# --- context per group

determine_context <- function(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA) {
  contexts <- rep("Other", times = length(None))
  
  contexts[SD > 0] <- "Seg. Dup."
  contexts[TRF > 0] <- "Tandem Repeat"
  contexts[SR > 0] <- "Simple Repeat"
  contexts[rDNA > 0] <- "rDNA"
  contexts[CenSat > 0] <- "CenSat"
  contexts[Telomere > 0] <- "Telomere"
  
  return(contexts)
}

context_list <- c("Other", "Simple Repeat", "Tandem Repeat", "Seg. Dup.", "CenSat", "Telomere", "rDNA")

grouped_variants <- read_rds("results/real/grouped_variants_real_500_v2.rds")
annotation <- read_tsv(annotation_real, col_names = TRUE)
grouped_variants <- grouped_variants %>%
  inner_join(annotation, by = "VARIANT")

grouped_variants <- grouped_variants %>%
  filter(ORIGIN == "real - 2025", COVERAGE == 32.51, SPANNING == "complete", ALN_TYPE == "global") %>%
  mutate(CONTEXT = determine_context(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA))

group_counts <- grouped_variants %>%
  group_by(GROUP, TECHNOLOGY, METHOD) %>%
  summarize(N = n(), .groups = "keep")

p_group_context <- grouped_variants %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(GROUP, CONTEXT, TECHNOLOGY, METHOD) %>%
  summarize(N1 = n(), .groups = "keep") %>%
  inner_join(group_counts, by = c("GROUP", "TECHNOLOGY", "METHOD")) %>%
  mutate(Fraction = N1 / N * 100) %>%
  ggplot(aes(x = CONTEXT, fill = METHOD, y = Fraction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(GROUP), cols = vars(TECHNOLOGY)) +
  custom_theme +
  scale_fill_manual(limits = method_order_real, values = my_colors) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
  )

# --- How much does each repeat context contribute to missed alleles?
grouped_variants_het <- read_rds("results/real/grouped_variants_real_500_v2.rds")
annotation <- read_tsv(annotation_real, col_names = TRUE)
grouped_variants_het <- grouped_variants_het %>%
  inner_join(annotation, by = "VARIANT")

grouped_variants_het <- grouped_variants_het %>%
  filter(ORIGIN == "real - 2025", COVERAGE == 32.51, SPANNING == "complete", ALN_TYPE == "global", !HOMOZYGOUS) %>%
  mutate(CONTEXT = determine_context(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA))


group_counts_het <- grouped_variants_het %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(GROUP, TECHNOLOGY, METHOD) %>%
  summarize(N = n(), .groups = "keep") %>%
  group_by(TECHNOLOGY, METHOD) %>%
  summarize(
    Missed = sum(N[GROUP == "E"]) + sum(N[GROUP == "F"]) + 2*sum(N[GROUP == "G"]),
    FP = sum(N[GROUP == "F"]),
    TP = 2*sum(N[GROUP == "D"]) + sum(N[GROUP == "E"]) + sum(N[GROUP == "F"]),
    .groups = "keep"
  )

p_group_context_het <- grouped_variants_het %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(GROUP, CONTEXT, TECHNOLOGY, METHOD) %>%
  summarize(N = n(), .groups = "keep") %>%
  group_by(TECHNOLOGY, METHOD, CONTEXT) %>%
  summarize(
    MissedInContext = sum(N[GROUP == "E"]) + sum(N[GROUP == "F"]) + 2*sum(N[GROUP == "G"]), .groups = "keep"
  ) %>%
  inner_join(group_counts, by = c("TECHNOLOGY", "METHOD")) %>%
  mutate(Fraction = MissedInContext / Missed * 100) %>%
  ggplot(aes(x = CONTEXT, fill = CONTEXT, y = Fraction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(METHOD), cols = vars(TECHNOLOGY)) +
  custom_theme +
  scale_fill_manual(values = my_colors) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
  )

grouped_variants_het %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(GROUP, CONTEXT, TECHNOLOGY, METHOD) %>%
  summarize(N = n(), .groups = "keep") %>%
  group_by(TECHNOLOGY, METHOD, CONTEXT) %>%
  summarize(
    MissedInContext = sum(N[GROUP == "E"]) + sum(N[GROUP == "F"]) + 2*sum(N[GROUP == "G"]), .groups = "keep"
  ) %>%
  inner_join(group_counts, by = c("TECHNOLOGY", "METHOD")) %>%
  mutate(Fraction = MissedInContext / Missed * 100) %>%
  write_tsv(paste0(numbers_dir, "/missed_fraction_per_context.tsv"))

pdf(paste0(output_dir, "/group_counts_context.pdf"), width = 20, height = 25)
plot(p_group_context)
plot(p_group_context_het)
dev.off()

# ---

determine_context <- function(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA) {
  contexts <- rep("Other", times = length(None))
  
  contexts[SD > 0] <- "Seg. Dup."
  contexts[TRF > 0] <- "Tandem Repeat"
  contexts[SR > 0] <- "Simple Repeat"
  contexts[rDNA > 0] <- "rDNA"
  contexts[CenSat > 0] <- "CenSat"
  contexts[Telomere > 0] <- "Telomere"
  
  return(contexts)
}

context_list <- c("Other", "Simple Repeat", "Tandem Repeat", "Seg. Dup.", "CenSat", "Telomere", "rDNA")

grouped_variants <- read_rds("results/real/grouped_variants_real_500_v2.rds") %>%
  bind_rows(read_rds("results/real/grouped_variants_real_1000_v2.rds"))
annotation <- read_tsv(annotation_real, col_names = TRUE)
grouped_variants <- grouped_variants %>%
  inner_join(annotation, by = "VARIANT")
grouped_variants <- grouped_variants %>%
  filter(COVERAGE == 32.51) %>%
  mutate(CONTEXT = determine_context(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA))

pr_values <- grouped_variants %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, CONTEXT, MARGIN, ORIGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, CONTEXT, MARGIN, ORIGIN) %>%
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
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal"))

pr_values$SPANNING[pr_values$SPANNING == "complete"] <- "spanning"
p_context_500_20 <- pr_values %>%
  filter(ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "spanning", MARGIN == 500, ORIGIN == "real - 2020", METHOD != "Multi POA v2") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(fill = factor(METHOD, levels = method_order_real), y = Value * 100, x = factor(CONTEXT, levels = rev(context_list)))) +
  facet_grid(rows = vars(Metric), cols = vars(interaction(TECHNOLOGY, SPANNING, sep = " - "))) +
  custom_theme +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "Repeat context")) +
  ggtitle("Flanking sequence size: 500 (2020)")

p_context_500_25 <- pr_values %>%
  filter(ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "spanning", MARGIN == 500, ORIGIN == "real - 2025", METHOD != "Multi POA v2") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(fill = factor(METHOD, levels = method_order_real), y = Value * 100, x = factor(CONTEXT, levels = rev(context_list)))) +
  facet_grid(rows = vars(Metric), cols = vars(interaction(TECHNOLOGY, SPANNING, sep = " - "))) +
  custom_theme +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "Repeat context")) +
  ggtitle("Flanking sequence size: 500 (2025)")

p_context_1000_20 <- pr_values %>%
  filter(ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "spanning", MARGIN == 1000, ORIGIN == "real - 2020", METHOD != "Multi POA v2") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(fill = factor(METHOD, levels = method_order_real), y = Value * 100, x = factor(CONTEXT, levels = rev(context_list)))) +
  facet_grid(rows = vars(Metric), cols = vars(interaction(TECHNOLOGY, SPANNING, sep = " - "))) +
  custom_theme +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "Repeat context")) +
  ggtitle("Flanking sequence size: 1000 (2020)")

p_context_1000_25 <- pr_values %>%
  filter(ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "spanning", MARGIN == 1000, ORIGIN == "real - 2025", METHOD != "Multi POA v2") %>%
  pivot_longer(cols = c("Precision", "Recall"), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(fill = factor(METHOD, levels = method_order_real), y = Value * 100, x = factor(CONTEXT, levels = rev(context_list)))) +
  facet_grid(rows = vars(Metric), cols = vars(interaction(TECHNOLOGY, SPANNING, sep = " - "))) +
  custom_theme +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(fill = guide_legend(title = "Repeat context")) +
  ggtitle("Flanking sequence size: 1000 (2025)")

pdf(paste0(output_dir, "/repeat_context.pdf"), width = 20, height = 8)
plot(p_context_500_20)
plot(p_context_500_25)
plot(p_context_1000_20)
plot(p_context_1000_25)
dev.off() 

grouped_variants %>%
  group_by(CONTEXT, MARGIN, ORIGIN) %>%
  summarize(N = length(unique(VARIANT)), .groups = "keep") %>%
  write_tsv(paste0(numbers_dir, "/variant_context_numbers.tsv"))

pr_values %>%
  filter(ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "complete") %>%
  write_tsv(paste0(numbers_dir, "/variant_context_scores.tsv"))

# try error rates per context
true_positives <- read_rds("results/real/true_positives_real_v2.rds") %>%
  filter(COVERAGE == 32.51, ALN_TYPE == "local" & SPANNING == "partial" | ALN_TYPE == "global" & SPANNING == "complete") %>%
  inner_join(annotation, by = "VARIANT") %>%
  mutate(CONTEXT = determine_context(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA))
# true_positives$SPANNING[true_positives$SPANNING == "complete"] <- "spanning"

p_rep_er_20 <- true_positives %>%
  filter(ORIGIN == "real - 2020", METHOD != "Multi POA v2") %>%
  group_by(TECHNOLOGY, SPANNING, CONTEXT, MARGIN, METHOD) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    # seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot() +
  facet_grid(rows = vars(interaction(TECHNOLOGY, SPANNING, sep = " - ")), cols = vars(MARGIN)) +
  geom_bar(
    aes(x = factor(CONTEXT, levels = rev(context_list)), fill = factor(METHOD, levels = method_order_real), y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  # geom_errorbar(
  #   aes(ymin = AvgER - seER, ymax = AvgER + seER, col = factor(METHOD, levels = method_order_real)),
  #   position = "dodge"
  # ) +
  guides(
    fill = guide_legend(title = "Method"),
    col = guide_none()
  ) +
  custom_theme + 
  ylab("Median Error Rate [%]") +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  scale_color_manual(values = my_colors, limits = method_order_real) +
  theme(
    axis.title.x = element_blank()
  ) +
  theme(
    panel.grid.major.x = element_blank()
  ) +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9)
  ) +
  ggtitle("Error Rate by sequence context - 2020 data")

p_rep_er_25 <- true_positives %>%
  filter(ORIGIN == "real - 2025", METHOD != "Multi POA v2") %>%
  group_by(TECHNOLOGY, SPANNING, CONTEXT, MARGIN, METHOD) |>
  summarize(
    AvgER = mean(ErrorRate * 100),
    # seER = sd(ErrorRate * 100) / sqrt(n()),
    .groups = "keep"
  ) |>
  ggplot() +
  facet_grid(rows = vars(interaction(TECHNOLOGY, SPANNING, sep = " - ")), cols = vars(MARGIN)) +
  geom_bar(
    aes(x = factor(CONTEXT, levels = rev(context_list)), fill = factor(METHOD, levels = method_order_real), y = AvgER), 
    position = "dodge", stat = "identity"
  ) +
  # geom_errorbar(
  #   aes(ymin = AvgER - seER, ymax = AvgER + seER, col = factor(METHOD, levels = method_order_real)),
  #   position = "dodge"
  # ) +
  guides(
    fill = guide_legend(title = "Method"),
    col = guide_none()
  ) +
  custom_theme + 
  ylab("Median Error Rate [%]") +
  scale_fill_manual(values = my_colors, limits = method_order_real) +
  scale_color_manual(values = my_colors, limits = method_order_real) +
  theme(
    axis.title.x = element_blank()
  ) +
  theme(
    panel.grid.major.x = element_blank()
  ) +
  theme(
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9)
  ) +
  ggtitle("Error Rate by sequence context - 2025 data")


pdf(paste0(output_dir, "/repeat_context_error_rate.pdf"), width = 20, height = 13)
plot(p_rep_er_20)
plot(p_rep_er_25)
dev.off()

# --- F1 score by method and SV type for simulated data

group_counts_sim <- read_rds(paste0(results_dir,"/grouped_variants.rds"))

group_counts_sim <- group_counts_sim %>%
  filter(METHOD != "Multi POA v2", ALN_TYPE == "global", SPANNING == "complete", COVERAGE == 30, ALN_PARAMS == "1_-1_-1_-1", ACCURACY == 0.99)

group_counts_sim <- group_counts_sim %>% mutate(SVTYPE = sapply(strsplit(VARIANT, split = "_"), function(x) {return(paste(x[1:(length(x) - 2)], collapse = "_"))})) 

pr_sim <- group_counts_sim %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, COVERAGE, SVTYPE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, COVERAGE, SVTYPE) %>%
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
  mutate(ORIGIN = "simulation")

p_f1_sim <- pr_sim %>%
  ggplot(aes(x = METHOD, y = SVTYPE)) +
  # geom_tile(fill = "white",colour = "grey", alpha = 0.2) +
  geom_point(aes(col = F1, size = log2(F1))) +
  facet_grid(cols = vars(TECHNOLOGY)) +
  scale_size_continuous(range = c(1, 20), limits = c(-1,0)) +
  scale_color_gradient2(limits = c(0,1), midpoint = 0.5, low = "black", mid = "blue", high = "green") +
  guides(
    size = guide_none()
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(limits = c("True phase", "A priori", "abPOA", "Multi POA", "RapidFuzz", "k-mer clust"))

p_f1_sim_v2 <- pr_sim %>%
  group_by(SVTYPE, TECHNOLOGY, COVERAGE, ORIGIN) %>%
  summarize(AvgF1 = mean(F1), se = sqrt(sd(F1)/6), .groups = "keep") %>%
  ggplot(aes(x = TECHNOLOGY)) +
  geom_bar(aes(fill = SVTYPE, y = AvgF1), stat = "identity", position = "dodge", col = "white") +
  geom_errorbar(aes(ymin = AvgF1 - se, ymax = AvgF1 + se, col = SVTYPE), position = "dodge") +
  custom_theme +
  ylab("Average F1 score") +
  theme(axis.title.x = element_blank())

pdf(paste0(output_dir, "/f1_sim_svtype.pdf"), width = 24, height = 10)
plot(p_f1_sim)
dev.off()

pdf(paste0(output_dir, "/f1_sim_svtype_avg.pdf"), width = 16, height = 9)
plot(p_f1_sim_v2)
dev.off()

pr_sim %>%
  write_tsv(paste0(numbers_dir, "/f1_sim_type.tsv"))

# --- plot precision and recall for each method across different alignment parameters in the simulation
pr_sim <- read_rds(paste0(results_dir, "/pr_values.rds"))
p_f1_aln_v1 <- pr_sim %>% 
  filter((ALN_TYPE == "global" & SPANNING == "complete") | (ALN_TYPE == "local" & SPANNING == "partial"), COVERAGE == 30, ACCURACY == 0.99, METHOD != "Multi POA v2") %>%
  mutate(SPANNING = ifelse(SPANNING == "complete", "spanning", "partial")) %>%
  pivot_longer(cols = c("Recall", "Precision"), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Value, y = interaction(TECHNOLOGY, SPANNING, ALN_PARAMS, sep = " - "), fill = METHOD)) +
  geom_bar(position = "dodge", stat = "identity", col = "black") +
  facet_grid(cols = vars(Metric)) +
  custom_theme +
  theme(
    axis.title.y = element_blank()
  )
  
p_f1_aln_v2 <- pr_sim %>% 
  filter((ALN_TYPE == "global" & SPANNING == "complete") | (ALN_TYPE == "local" & SPANNING == "partial"), COVERAGE == 30, ACCURACY == 0.99, METHOD != "Multi POA v2") %>%
  mutate(SPANNING = ifelse(SPANNING == "complete", "spanning", "partial")) %>%
  group_by(TECHNOLOGY, SPANNING, ALN_PARAMS) %>%
  summarize(AvgF1 = mean(F1), seF1 = sqrt(sd(F1)/6), .groups = "keep") %>%
  ggplot(aes(x = AvgF1, y = interaction(TECHNOLOGY, SPANNING, ALN_PARAMS, sep = " - "), fill = ALN_PARAMS)) +
  geom_bar(position = "dodge", stat = "identity", col = "black") +
  geom_errorbar(aes(xmin = AvgF1 - seF1, xmax = AvgF1 + seF1, col = ALN_PARAMS), position = "dodge") +
  custom_theme +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("F1 score")

pdf(paste0(output_dir, "/sim_f1_by_aln_params_detailed.pdf"), width = 16, height = 16)  
plot(p_f1_aln_v1)
dev.off()

pdf(paste0(output_dir, "/sim_f1_by_aln_params.pdf"), width = 16, height = 9)  
plot(p_f1_aln_v2)
dev.off()


