source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

### load GIAB simulated data
variant_regions_giab <- read_tsv(variant_file_giab, col_names = FALSE)
colnames(variant_regions_giab) <- c("VARIANT", "CHR", "TYPE", "GT", "REGION")
variant_regions_giab <- variant_regions_giab %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = zygosity_from_gt(GT)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM")

no_overlap_giab <- read_tsv(giab_overlap_whitelist, col_names = "VARIANT")
no_low_cov <- read_tsv(variant_coverage_whitelist_old)
variant_regions_giab <- variant_regions_giab %>%
  inner_join(no_overlap_giab, by = "VARIANT") %>% 
  inner_join(no_low_cov, by = "VARIANT")

variant_regions_giab %>% group_by(TYPE) %>% summarize(N = n())

results_giab <- tibble()
for (technology in c("ont", "pacbio"))
{
  cov_dirs <- list.dirs(paste0(results_dir_giab, "/", technology), recursive = FALSE, full.names = TRUE)
  for (cov_dir in cov_dirs) {
    cov <- as.numeric(tail(strsplit(cov_dir, split = "/")[[1]], n = 1))
    # extract coverage
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
      
      
      method <- gsub(info, pattern = "evaluation_", replacement = "")
      cat("\nMethod: ", method, "\n\n")
      
      results_giab <- results_giab |> rbind(
        read_tsv(ef, col_names = TRUE) |>
          mutate(VARIANT = factor(VARIANT, levels = variant_regions_giab$VARIANT)) |>
          complete(VARIANT) |>
          mutate(
            COVERAGE = cov,
            TECHNOLOGY = technology,
            SPANNING = spanning,
            ALN_TYPE = aln_type,
            METHOD = method
          )
      )
    }
  }
}

# no longer needed, but created at some point during testing
results_giab <- results_giab %>%
  filter(
    METHOD != "abPOA_conservative", 
    METHOD != "fuzzy_v2", 
    METHOD != "fuzzy_hclust", 
    METHOD != "ITER_POA_ONT",
    METHOD != "fuzzy_v3",
    METHOD != "ITER_POA_INTERMEDIATE"#,
    # METHOD != "ITER_POA_PERMISSIVE"
  )

# calculate alignment match rate and join TYPE and ZYGOSITY
results_giab <- results_giab %>% 
  mutate(MATCHED_BASES = as.numeric(MATCHED_BASES), SCORE = as.numeric(SCORE), ALLELE_MATCH = as.numeric(ALLELE_MATCH), CONSENSUS_LENGTH = as.numeric(CONSENSUS_LENGTH), CONSENSUS_IDX = as.numeric(CONSENSUS_IDX)) %>%
  mutate(MATCH_RATE = MATCHED_BASES / CONSENSUS_LENGTH) %>%
  inner_join(variant_regions_giab, by = c("VARIANT")) 

results_giab <- results_giab %>%
  mutate(TECHNOLOGY = toupper(TECHNOLOGY))
results_giab[which(results_giab$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"

results_giab <- results_giab %>%
  mutate(METHOD = rename_methods(METHOD)) %>%
  mutate(MARGIN = 500)

## same preprocessing as for real data
nvar <- results_giab %>% pull(VARIANT) %>% unique() %>% length

p1 <- results_giab %>% 
  filter(ALLELE_MATCH >= 0, MATCH_RATE < 0.95, COVERAGE == 15, ALN_TYPE == "global") %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>%
  ggplot(aes(x = TECHNOLOGY, fill = factor(METHOD, levels = method_order), y = N / nvar * 100)) +
  geom_bar(position = "dodge", stat = "identity") +
  custom_theme +
  ylab("Variant Loci [%]") +
  scale_fill_manual(limits = method_order, values = my_colors) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )

p2 <- results_giab %>% 
  filter(ALLELE_MATCH >= 0, MATCH_RATE < 0.95, COVERAGE == 15, ALN_TYPE == "global") %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>%
  ggplot(aes(x = TECHNOLOGY, fill = factor(METHOD, levels = method_order), y = N)) +
  geom_bar(position = "dodge", stat = "identity") +
  custom_theme +
  ylab("Variant Loci") +
  scale_fill_manual(limits = method_order, values = my_colors)+
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
pdf(paste0(output_dir, "/high_error_loci_giab.pdf"), width = 16, height = 9)
plot(p1)
plot(p2)
dev.off()

results_giab %>% 
  filter(COVERAGE == 15, ALN_TYPE == "global",ALLELE_MATCH >= 0, MATCH_RATE < 0.95) %>%
  filter(METHOD != "Multi POA v2") %>%
  group_by(METHOD, TECHNOLOGY, VARIANT) %>%
  summarize(N_distinct = n(), .groups = "keep") %>%
  group_by(METHOD, TECHNOLOGY) %>%
  summarize(N = n(), .groups = "keep") %>% write_rds(paste0(numbers_dir, "/high_error_loci_numbers_giab.tsv"))


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

context_assignments <- results_giab %>%
  filter(ALLELE_MATCH >= 0, MATCH_RATE < 0.95) %>%
  inner_join(read_tsv(annotation_real), by = "VARIANT") %>%
  filter(COVERAGE == 15, SPANNING == "complete", ALN_TYPE == "global") %>%
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
pdf(paste0(output_dir, "/high_error_loci_giab_context.pdf"), width = 16, height = 9)
plot(p3)
dev.off()



results_giab$ALLELE_MATCH[results_giab$MATCH_RATE < 0.95] <- -1
results_giab %>% write_rds(paste0(results_dir_giab, "/results_giab_full.rds"))


## Group variants based on allele assignments
grouped_variants_giab <- results_giab %>%
  group_by(
    VARIANT, METHOD, COVERAGE, TECHNOLOGY, HOMOZYGOUS, ALN_TYPE, SPANNING, MARGIN
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
  mutate(GROUP = determine_group(HOMOZYGOUS, N, N_unique, NoSeq, ALL_UNCOVERED, ONE_UNCOVERED, ALL_UNMATCHED, ONE_UNMATCHED)) %>%
  mutate(COVERAGE = 2 * COVERAGE)
grouped_variants_giab %>% write_rds(paste0(results_dir_giab, "/grouped_variants_giab_full.rds"))

pr_values_giab <- grouped_variants_giab %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN) %>%
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
  filter(!(METHOD == "abPOA" & ALN_TYPE == "semiglobal")) %>%
  mutate(ORIGIN = "giab_sim")
pr_values_giab %>% write_rds(paste0(results_dir_giab, "/pr_values_giab_full.rds"))

true_positives <- results_giab %>%
  filter(!is.na(MATCH_RATE), !is.na(BP_COVERED), BP_COVERED == 1, ALLELE_MATCH != -1) %>%
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, MARGIN) |>
  summarize(
    ErrorRate = min(1 - MATCH_RATE), 
    .groups = "keep"
  ) |>
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, MARGIN) |>
  summarize(
    ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
    .groups = "keep"
  )
true_positives %>% write_rds(paste0(results_dir_giab, "/true_positives_full.rds"))

rm(list = ls())
gc()
