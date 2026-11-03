source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")

### process data from main simulation 

# load variant info, esp. zygosity and type
variant_regions <- read_tsv(variant_file, col_names = FALSE)
colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "ZYGOSITY", "REGION")
variant_regions <- variant_regions %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = toupper(ZYGOSITY)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM")

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
          spanning <- ifelse(grepl(ef, pattern = "complete") || (aln_type == "global"), "complete", "partial")
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
          info <- gsub(info, pattern = "_local", replacement = "")
          info <- gsub(info, pattern = "_partial", replacement = "")
          info <- gsub(info, pattern = "_complete", replacement = "")
          
          
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
    METHOD != "ITER_POA_INTERMEDIATE"#,
    #METHOD != "ITER_POA_PERMISSIVE"
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
  mutate(MATCHED_BASES = as.numeric(MATCHED_BASES), SCORE = as.numeric(SCORE), ALLELE_MATCH = as.numeric(ALLELE_MATCH), CONSENSUS_LENGTH = as.numeric(CONSENSUS_LENGTH), CONSENSUS_IDX = as.numeric(CONSENSUS_IDX)) %>%
  mutate(MATCH_RATE = MATCHED_BASES / CONSENSUS_LENGTH) %>%
  inner_join(variant_regions, by = c("VARIANT")) %>%
  mutate(COMPHET = grepl(TYPE, pattern = "COMPHET"))


results <- results %>%
  mutate(TECHNOLOGY = toupper(TECHNOLOGY))
results[which(results$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"

results <- results %>%
  mutate(METHOD = rename_methods(METHOD))

results$ALLELE_MATCH[results$MATCH_RATE < 0.95] <- -1
write_rds(results, paste0(results_dir, "/results.rds"))

true_positives <- results %>%
  filter(!is.na(MATCH_RATE), !is.na(BP_COVERED), BP_COVERED == 1) %>%
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, ALLELE_MATCH, ACCURACY, ALN_PARAMS, VARSIZE) |>
  summarize(
    ErrorRate = min(1 - MATCH_RATE), 
    .groups = "keep"
  ) |>
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, ACCURACY, ALN_PARAMS, VARSIZE) |>
  summarize(
    ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
    .groups = "keep"
  )
true_positives %>% write_rds(paste0(results_dir, "/true_positives.rds"))
rm(true_psoitives)
gc()

grouped_variants <- results %>%
  group_by(
    VARIANT, METHOD, ALN_PARAMS, COVERAGE, TECHNOLOGY, ACCURACY, VARSIZE, HOMOZYGOUS, ALN_TYPE, SPANNING
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
  mutate(COVERAGE = 2*COVERAGE)
grouped_variants %>% write_rds(paste0(results_dir, "/grouped_variants.rds"))
rm(results)
gc()

grouped_variants <- read_rds(paste0(results_dir, "/grouped_variants.rds"))
pr_values <- grouped_variants %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, ALN_PARAMS, ACCURACY) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, ALN_PARAMS, ACCURACY) %>%
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
  mutate(ORIGIN = "simulation")

pr_values %>% write_rds(paste0(results_dir, "/pr_values.rds"))



