source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

### load unertainty simulation data
results_sim <- tibble()
for (margin in c(100, 500, 1000))
{
  for (jitter in c(0, 50, 200))
  {
    # load variant info, esp. zygosity and type
    # Variant region file depends on the result i am looking at
    variant_file <- paste0("data/simulation/variant_regions_", margin, "_", jitter, ".tsv")
    variant_regions <- read_tsv(variant_file, col_names = FALSE)
    colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "ZYGOSITY", "REGION")
    variant_regions <- variant_regions %>%
      select(-REGION, -CHR) %>%
      mutate(ZYGOSITY = toupper(ZYGOSITY)) %>%
      mutate(HOMOZYGOUS = ZYGOSITY == "HOM")
    
    for (technology in c("ont", "pacbio"))
    {
      eval_files <- list.files(
        paste0(results_dir_sim, technology, "/", margin, "/", jitter),
        full.names = TRUE,
        pattern = "evaluation_.*.tsv"
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
        
        temp_results <- read_tsv(ef, col_names = TRUE) |>
          mutate(VARIANT = factor(VARIANT, levels = variant_regions$VARIANT)) |>
          complete(VARIANT) |>
          mutate(
            COVERAGE = 15,
            TECHNOLOGY = technology,
            SPANNING = spanning,
            ALN_TYPE = aln_type,
            METHOD = method,
            ACCURACY = 0.99,
            MARGIN = margin,
            JITTER = jitter
          )
        temp_results <- temp_results %>%
          filter(
            METHOD != "abPOA_conservative", 
            METHOD != "fuzzy_v2", 
            METHOD != "fuzzy_hclust", 
            METHOD != "ITER_POA_ONT",
            METHOD != "fuzzy_v3",
            METHOD != "ITER_POA_INTERMEDIATE"#,
            # METHOD != "ITER_POA_PERMISSIVE"
          )
        temp_results <- temp_results %>% mutate(
          VARSIZE = sapply(
            as.character(VARIANT),
            function(x) {
              return(as.numeric(strsplit(x, split = "_")[[1]][length(strsplit(x, split = "_")[[1]])]))
            }
          )
        ) %>%
          mutate(MATCHED_BASES = as.numeric(MATCHED_BASES), SCORE = as.numeric(SCORE), ALLELE_MATCH = as.numeric(ALLELE_MATCH), CONSENSUS_LENGTH = as.numeric(CONSENSUS_LENGTH), CONSENSUS_IDX = as.numeric(CONSENSUS_IDX)) %>%
          mutate(MATCH_RATE = MATCHED_BASES / CONSENSUS_LENGTH) %>%
          inner_join(variant_regions, by = c("VARIANT")) 
        
        temp_results <- temp_results %>%
          mutate(TECHNOLOGY = toupper(TECHNOLOGY))
        temp_results[which(temp_results$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"
        results_sim <- results_sim %>% rbind(temp_results)
      }
    }
  }
}
results_sim <- results_sim %>% mutate(METHOD = rename_methods(METHOD))
results_sim$ALLELE_MATCH[results_sim$MATCH_RATE < 0.95] <- -1
results_sim %>% write_rds(paste0(results_dir_sim, "/results_uncertainty.rds"))

grouped_variants_sim <- results_sim %>%
  group_by(
    VARIANT, METHOD, COVERAGE, TECHNOLOGY, HOMOZYGOUS, ALN_TYPE, SPANNING, MARGIN, JITTER
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
grouped_variants_sim %>% write_rds(paste0(results_dir_sim, "/grouped_variants_uncertainty.rds"))

pr_values_sim <- grouped_variants_sim %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, JITTER) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, MARGIN, JITTER) %>%
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
  mutate(ORIGIN = "uncertainty_sim") %>%
  mutate(COVERAGE = 2*COVERAGE)

pr_values_sim %>% write_rds(paste0(results_dir_sim, "/pr_values_uncertainty.rds"))

true_positives <- results_sim %>%
  filter(!is.na(MATCH_RATE), !is.na(BP_COVERED), BP_COVERED == 1) %>%
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, ALLELE_MATCH, MARGIN, JITTER) |>
  summarize(
    ErrorRate = min(1 - MATCH_RATE), 
    .groups = "keep"
  ) |>
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, MARGIN, JITTER) |>
  summarize(
    ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
    .groups = "keep"
  )
true_positives %>% write_rds(paste0(results_dir_sim, "/true_positives.rds"))

