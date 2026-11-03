source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

### load variant region info, stats, annotations and remove low-coverage regions and overlapping regions ###

# paths
variant_file_real <- "data/references/HG002_assembly/variant_regions_500.txt"
results_dir_real <- "results/real/"

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

variant_regions <- variant_regions %>% 
  inner_join(no_low_cov, by = "VARIANT")

# reduce to regions that do not overlap any other regions
no_overlap_500 <- read_tsv(variant_overlap_whitelist_500, col_names = c("VARIANT"))
variant_regions_500 <- variant_regions %>%
  inner_join(no_overlap_500, by = "VARIANT")
variant_regions_500 %>% group_by(HOMOZYGOUS) %>% summarize(N = n())

no_overlap_1000 <- read_tsv(variant_overlap_whitelist_1000, col_names = c("VARIANT"))
variant_regions_1000 <- variant_regions %>%
  inner_join(no_overlap_1000, by = "VARIANT")
variant_regions_1000 %>% group_by(HOMOZYGOUS) %>% summarize(N = n())

### load results from real dataset 
# "ont_old", "pacbio_old"

results <- tibble()
for (technology in c("ont", "pacbio"))
{
  cov_dirs <- list.dirs(paste0(results_dir_real, "/", technology), recursive = FALSE, full.names = TRUE)
  for (cov_dir in cov_dirs) {
    cov <- as.numeric(tail(strsplit(cov_dir, split = "/")[[1]], n = 1))
    
    if (is.na(cov)) {
      next
    }
    eval_files <- list.files(
      cov_dir,
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
      spanning <- ifelse(grepl(ef, pattern = "complete") || (aln_type == "global"), "complete", "partial")
      
      info <- strsplit(ef, split = "/")[[1]]
      info <- gsub(info[length(info)], pattern = ".tsv", replacement = "")
      info <- gsub(info, pattern = "_semiglobal", replacement = "")
      info <- gsub(info, pattern = "_global", replacement = "")
      info <- gsub(info, pattern = "_complete", replacement = "")
      info <- gsub(info, pattern = "_partial", replacement = "")
      info <- gsub(info, pattern = "_local", replacement = "")
      
      info <- strsplit(info, split = "_")[[1]]
      info <- info[2:length(info)]
      method <- paste0(info[1:length(info) - 1], collapse = "_")
      margin <- as.integer(info[length(info)])
      cat("\nMethod: ", method, "\tMargin: ", margin, "\n\n")
      
      results <- results |> rbind(
        read_tsv(ef, col_names = TRUE) |>
          mutate(VARIANT = factor(VARIANT, levels = unique(variant_regions$VARIANT))) |>
          complete(VARIANT) |>
          mutate(
            COVERAGE = cov,
            TECHNOLOGY = technology,
            SPANNING = spanning,
            ALN_TYPE = aln_type,
            METHOD = method,
            MARGIN = margin,
            ORIGIN = "real - 2025"
          )
      )
    }
  }
}

for (technology in c("ont", "pacbio"))
{
  cov_dirs <- list.dirs(paste0(results_dir_real, "/", technology, "_old"), recursive = FALSE, full.names = TRUE)
  for (cov_dir in cov_dirs) {
    cov <- as.numeric(tail(strsplit(cov_dir, split = "/")[[1]], n = 1))
    
    if (is.na(cov)) {
      next
    }
    eval_files <- list.files(
      cov_dir,
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
      spanning <- ifelse(grepl(ef, pattern = "complete") || (aln_type == "global"), "complete", "partial")
      
      info <- strsplit(ef, split = "/")[[1]]
      info <- gsub(info[length(info)], pattern = ".tsv", replacement = "")
      info <- gsub(info, pattern = "_semiglobal", replacement = "")
      info <- gsub(info, pattern = "_global", replacement = "")
      info <- gsub(info, pattern = "_complete", replacement = "")
      info <- gsub(info, pattern = "_partial", replacement = "")
      info <- gsub(info, pattern = "_local", replacement = "")
      
      info <- strsplit(info, split = "_")[[1]]
      info <- info[2:length(info)]
      method <- paste0(info[1:length(info) - 1], collapse = "_")
      margin <- as.integer(info[length(info)])
      cat("\nMethod: ", method, "\tMargin: ", margin, "\n\n")
      
      results <- results |> rbind(
        read_tsv(ef, col_names = TRUE) |>
          mutate(VARIANT = factor(VARIANT, levels = unique(variant_regions$VARIANT))) |>
          complete(VARIANT) |>
          mutate(
            COVERAGE = cov,
            TECHNOLOGY = technology,
            SPANNING = spanning,
            ALN_TYPE = aln_type,
            METHOD = method,
            MARGIN = margin,
            ORIGIN = "real - 2020"
          )
      )
    }
  }
}

# results <- results %>%
#   filter(METHOD != "ITER_POA_PERMISSIVE")
results <- results %>%
  mutate(METHOD = rename_methods(METHOD))
results <- results %>% mutate(TECHNOLOGY = toupper(TECHNOLOGY))

# results %>% pull(ORIGIN) %>% table
results[which(results$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"
# results[which(results$TECHNOLOGY == "PACBIO_OLD"), "TECHNOLOGY"] <- "PacBio HiFi - 2020"
# results[which(results$TECHNOLOGY == "ONT_OLD"), "TECHNOLOGY"] <- "ONT - 2020"


results <- results %>% 
  mutate(
    MATCHED_BASES = as.numeric(MATCHED_BASES), 
    SCORE = as.numeric(SCORE), 
    ALLELE_MATCH = as.numeric(ALLELE_MATCH), 
    CONSENSUS_LENGTH = as.numeric(CONSENSUS_LENGTH), 
    CONSENSUS_IDX = as.numeric(CONSENSUS_IDX)
  ) %>%
  mutate(MATCH_RATE = MATCHED_BASES / CONSENSUS_LENGTH) %>%
  inner_join(variant_regions, by = c("VARIANT"))



nvar <- results %>% filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT")) %>% pull(VARIANT) %>% unique() %>% length

p1 <- results %>% 
  filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT")) %>%
  filter(MATCH_RATE < 0.95, ALLELE_MATCH >= 0) %>%
  filter(COVERAGE == 32.51, ALN_TYPE == "global", ORIGIN == "real - 2025") %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>%
  ggplot(aes(x = TECHNOLOGY, fill = factor(METHOD, levels = method_order_real), y = N / nvar * 100)) +
  geom_bar(position = "dodge", stat = "identity") +
  custom_theme +
  ylab("Variant Loci [%]") +
  scale_fill_manual(limits = method_order_real, values = my_colors) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

p2 <- results %>% 
  filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT")) %>%
  filter(MATCH_RATE < 0.95, ALLELE_MATCH >= 0) %>%
  filter(COVERAGE == 32.51, ALN_TYPE == "global", ORIGIN == "real - 2025") %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>%
  ggplot(aes(x = TECHNOLOGY, fill = factor(METHOD, levels = method_order_real), y = N)) +
  geom_bar(position = "dodge", stat = "identity") +
  custom_theme +
  ylab("Variant Loci") +
  scale_fill_manual(limits = method_order_real, values = my_colors)+
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
pdf(paste0(output_dir, "/high_error_loci_real.pdf"), width = 16, height = 9)
plot(p1)
plot(p2)
dev.off()

results %>% 
  filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT")) %>%
  filter(MATCH_RATE < 0.95, ALLELE_MATCH >= 0) %>%
  filter(COVERAGE == 32.51, ALN_TYPE == "global", ORIGIN == "real - 2025") %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>% write_rds(paste0(numbers_dir, "/high_error_loci_numbers_real.tsv"))


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

context_assignments <- results %>%
  filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT")) %>%
  filter(MATCH_RATE < 0.95, ALLELE_MATCH >= 0) %>%
  filter(ORIGIN == "real - 2025", COVERAGE == 32.51, SPANNING == "complete", ALN_TYPE == "global") %>%
  mutate(CONTEXT = determine_context(None, RM, SR, Alu, SD, CenSat, Telomere, TRF, rDNA))

p3 <- context_assignments %>%
  filter(METHOD != "Multi POA v2")%>%
  group_by(METHOD, TECHNOLOGY, VARIANT, CONTEXT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY, CONTEXT) %>%
  summarize(N = n(), .groups = "keep") %>%
  ggplot(aes(x = METHOD, fill = CONTEXT, y = N)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(rows = vars(TECHNOLOGY)) +
  custom_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9)
  ) +
  ylab("# consensus sequences")
pdf(paste0(output_dir, "/high_error_loci_real_context.pdf"), width = 16, height = 9)
plot(p3)
dev.off()





results$ALLELE_MATCH[results$MATCH_RATE < 0.95] <- -1


# results %>% write_rds(paste0(results_dir_real, "/real_results_all.rds"))
# results <- read_rds(paste0(results_dir_real, "/real_results_all.rds")) %>% filter(MARGIN == 500)
# results_all_500 <- results %>% filter(MARGIN == 500) %>% inner_join(no_overlap_500, by = "VARIANT")
# results_all_500 %>% write_rds(paste0(results_dir_real, "/real_results_all_500.rds"))
# rm(results_all_500)
# results_all_1000 <- results %>% filter(MARGIN == 1000) %>% inner_join(no_overlap_1000, by = "VARIANT")
# results_all_1000 %>% write_rds(paste0(results_dir_real, "/real_results_all_1000.rds"))
# 
# results_all_500 <- read_rds(paste0(results_dir_real, "/real_results_all_500.rds"))
# results_all_500 <- results_all_500 %>% inner_join(no_low_cov, by = "VARIANT")
# results_all_500 %>% write_rds(paste0(results_dir_real, "/real_results_500_v2.rds"))
# 
# results_all_1000 <- read_rds(paste0(results_dir_real, "/real_results_all_1000.rds"))
# results_all_1000 <- results_all_1000 %>% inner_join(no_low_cov, by = "VARIANT")
# results_all_1000 %>% write_rds(paste0(results_dir_real, "/real_results_1000_v2.rds"))

results_500 <- results %>%
  filter(MARGIN == 500) %>%
  inner_join(no_overlap_500, by = c("VARIANT"))
results_1000 <- results %>%
  filter(MARGIN == 1000) %>%
  inner_join(no_overlap_1000, by = "VARIANT")
rm(results)
gc()

results_500 %>% write_rds(paste0(results_dir_real, "/real_results_500_v2.rds"))
results_1000 %>% write_rds(paste0(results_dir_real, "/real_results_1000_v2.rds"))

## Group variants based on allele assignments
grouped_variants_500 <- results %>%
  group_by(
    VARIANT, METHOD, COVERAGE, TECHNOLOGY, HOMOZYGOUS, ALN_TYPE, SPANNING, MARGIN, ORIGIN
  ) %>%
  summarize(
    N = n(),
    N_unique = length(unique(ALLELE_MATCH[ALLELE_MATCH >= 0 & !is.na(ALLELE_MATCH)])),
    NoSeq = n() == 1 & any(is.na(ALLELE_MATCH)),
    ONE_UNCOVERED = any(BP_COVERED == 0),
    ALL_UNCOVERED = all(BP_COVERED == 0),
    ONE_UNMATCHED = any(ALLELE_MATCH == -1),
    ALL_UNMATCHED = all(ALLELE_MATCH == -1),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(GROUP = determine_group(HOMOZYGOUS, N, N_unique, NoSeq, ALL_UNCOVERED, ONE_UNCOVERED, ALL_UNMATCHED, ONE_UNMATCHED))

grouped_variants_500 %>% write_rds(paste0(results_dir_real, "/grouped_variants_real_500_v2.rds"))

pr_values_500 <- grouped_variants_500 %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
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

pr_values_500 %>% write_rds(paste0(results_dir_real, "/pr_values_real_500_v2.rds"))
rm(results_500)
rm(grouped_variants_500)
gc()

#results_1000 <- read_rds(paste0(results_dir_real, "/real_results_1000.rds"))
grouped_variants_1000 <- results_1000 %>%
  group_by(
    VARIANT, METHOD, COVERAGE, TECHNOLOGY, HOMOZYGOUS, ALN_TYPE, SPANNING, MARGIN, ORIGIN
  ) %>%
  summarize(
    N = n(),
    N_unique = length(unique(ALLELE_MATCH[ALLELE_MATCH >= 0 & !is.na(ALLELE_MATCH)])),
    NoSeq = n() == 1 & any(is.na(ALLELE_MATCH)),
    ONE_UNCOVERED = any(BP_COVERED == 0),
    ALL_UNCOVERED = all(BP_COVERED == 0),
    ONE_UNMATCHED = any(ALLELE_MATCH == -1),
    ALL_UNMATCHED = all(ALLELE_MATCH == -1),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(GROUP = determine_group(HOMOZYGOUS, N, N_unique, NoSeq, ALL_UNCOVERED, ONE_UNCOVERED, ALL_UNMATCHED, ONE_UNMATCHED))
grouped_variants_1000 %>% write_rds(paste0(results_dir_real, "/grouped_variants_real_1000_v2.rds"))

pr_values_1000 <- grouped_variants_1000 %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, ORIGIN) %>%
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
pr_values_1000 %>% write_rds(paste0(results_dir_real, "/pr_values_real_1000_v2.rds"))

pr_values_real <- pr_values_500 %>%
  mutate(MARGIN = 500) %>%
  bind_rows(pr_values_1000 %>% mutate(MARGIN = 1000))
pr_values_real %>% write_rds(paste0(results_dir_real, "/pr_values_real_v2.rds"))

rm(results_1000)
rm(grouped_variants_1000)
gc()

results <- read_rds(paste0(results_dir_real, "/real_results_500_v2.rds")) %>%
  mutate(MARGIN = 500) %>%
  bind_rows(
    read_rds(paste0(results_dir_real, "/real_results_1000_v2.rds")) %>%
      mutate(MARGIN = 1000)
  )

true_positives <- results %>%
  filter(!is.na(MATCH_RATE), !is.na(BP_COVERED), BP_COVERED == 1, ALLELE_MATCH != -1) %>%
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, MARGIN, ORIGIN) |>
  summarize(
    ErrorRate = min(1 - MATCH_RATE), 
    .groups = "keep"
  ) |>
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, MARGIN, ORIGIN) |>
  summarize(
    ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
    .groups = "keep"
  )
true_positives %>% write_rds(paste0(results_dir_real, "/true_positives_v2.rds"))

rm(list = ls())
gc()
