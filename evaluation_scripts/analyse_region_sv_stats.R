setwd("~/Git/2026_02_07_lr-consensus/")

source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

# differential_regions <- read_tsv("results/real_numbers/technology_loci_differences_32.51_complete_local.tsv")

region_stats <- read_tsv("data/real/HG002/variant_region_stats_ont_32.51.tsv") %>%
  mutate(TECHNOLOGY = "ONT", ORIGIN = "real") %>%
  bind_rows(
    read_tsv("data/real/HG002/variant_region_stats_pacbio_32.51.tsv") %>% 
      mutate(TECHNOLOGY = "PacBio HiFi", ORIGIN = "real")
  ) %>%
  bind_rows(
    read_tsv("data/simulation/giab/variant_region_stats_ont_15.tsv") %>%
      mutate(TECHNOLOGY = "ONT", ORIGIN = "sim")
  ) %>%
  bind_rows(
    read_tsv("data/simulation/giab/variant_region_stats_pacbio_15.tsv") %>%
      mutate(TECHNOLOGY = "PacBio HiFi", ORIGIN = "sim")
  )
variant_regions <- read_tsv(variant_file_real, col_names = c("VARIANT", "CHR", "TYPE", "GT", "REGION"))



non_overlapping <- read_tsv(variant_overlap_whitelist_500, col_names = "VARIANT")

region_stats <- region_stats %>%
  inner_join(variant_regions, by = "VARIANT") %>%
  inner_join(non_overlapping, by = "VARIANT") %>%
  mutate(ZYG = zygosity_from_gt(GT)) %>%
  mutate(HOMOZYGOUS = ZYG == "HOM")


region_stats %>%
  filter(ORIGIN == "real") %>%
  mutate(NO_SIGNALS = N_BP < 3) %>%
  group_by(HOMOZYGOUS, TECHNOLOGY, NO_SIGNALS) %>%
  summarize(N = n())


region_stats %>%
  group_by(TECHNOLOGY, ORIGIN) %>%
  summarize(low_var = sum(N_BP < 3) / n())

region_stats %>%
  group_by(TECHNOLOGY, ORIGIN) %>%
  summarize(low_cov = sum(COVERAGE < 6) / n())

region_stats %>%
  filter(ORIGIN == "real") %>%
  filter(COVERAGE > 5, N_BP >= 3) %>%
  group_by(VARIANT) %>%
  summarize(ALL = n() == 2) %>%
  filter(ALL) %>%
  select(-ALL) %>% write_tsv(variant_coverage_whitelist)


region_stats %>%
  filter(ORIGIN == "real") %>%
  filter(COVERAGE <= 5, N_BP < 3) %>%
  group_by(TECHNOLOGY) %>%
  summarise(N = n())

no_low_cov <- region_stats %>%
  filter(ORIGIN == "real") %>%
  filter(COVERAGE > 5, N_BP >= 3) %>%
  group_by(VARIANT) %>%
  summarize(ALL = n() == 2) %>%
  filter(ALL) %>%
  select(-ALL) %>% select(VARIANT)
region_stats <- region_stats %>% inner_join(no_low_cov, by = "VARIANT")

region_stats %>% group_by(TECHNOLOGY, ORIGIN, TYPE) %>% summarize(N = n())

p1 <- region_stats %>% 
  ggplot(aes(
    x = SD_POS, 
    group = ifelse(grepl(VARIANT, pattern = "INS"), "Ins", "Del"), 
    fill = ifelse(grepl(VARIANT, pattern = "INS"), "Ins", "Del")
    )) +
  geom_histogram(position = "dodge", bins = 100) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(ORIGIN)) +
  guides(fill = guide_legend(title = "SV Type")) +
  theme_bw()+
  custom_theme +
  xlab("SD (position)")

p3 <- region_stats %>%
  mutate(SVTYPE = ifelse(grepl(VARIANT, pattern = "DEL"), "DEL", "INS")) %>%
  pivot_longer(cols = c(SD_DEL_SIZE, SD_INS_SIZE), values_to = "SIZE_SD", names_to = "CATEGORY") %>% 
  mutate(SD_CAT_TYPE = ifelse(grepl(CATEGORY, pattern = "DEL"), "DEL", "INS")) %>%
  filter(SVTYPE == SD_CAT_TYPE) %>%
  ggplot(aes(
    x = SIZE_SD,
    group = SVTYPE,
    fill = SVTYPE
  )) +
  geom_histogram(position = "dodge", bins = 100) +
  facet_grid(rows = vars(TECHNOLOGY), cols = vars(ORIGIN)) +
  theme_bw() +
  custom_theme +
  xlab("SD (size)")
  

pdf(paste0(output_dir, "/technology_sv_stats.pdf"), width = 12, height = 7)
plot(p1)
plot(p1 + scale_y_log10() + ylab("Log count"))
plot(p3)
plot(p3 + scale_y_log10() + ylab("Log count"))
dev.off()
