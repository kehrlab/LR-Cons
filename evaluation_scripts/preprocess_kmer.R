source("evaluation_scripts/definitions.R")
source("evaluation_scripts/functions.R")


# load variant info, esp. zygosity and type
variant_regions <- read_tsv(variant_file_kmeans, col_names = FALSE)
colnames(variant_regions) <- c("VARIANT", "CHR", "TYPE", "ZYGOSITY", "REGION")
variant_regions <- variant_regions %>%
  select(-REGION, -CHR) %>%
  mutate(ZYGOSITY = toupper(ZYGOSITY)) %>%
  mutate(HOMOZYGOUS = ZYGOSITY == "HOM")

# load evaluation results for all methods
results <- tibble()
for (technology in c("ont", "pacbio"))
{
  acc_dirs <- list.dirs(paste0(results_dir_kmeans, "/", technology), recursive = FALSE, full.names = TRUE)
  for (acc_dir in acc_dirs) {
    acc <- as.numeric(tail(strsplit(acc_dir, split = "/")[[1]], n = 1))
    cov_dirs <- list.dirs(acc_dir, recursive = FALSE, full.names = TRUE)
    for (cov_dir in cov_dirs)
    {
      cov <- as.numeric(tail(strsplit(cov_dir, split = "/")[[1]], n = 1))
      # extract coverage
      if (is.na(cov)) {
        next
      }
      eval_files <- list.files(
        paste0(cov_dir, "/1_-1_-1_-1/"),
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
        
        info <- gsub(info, pattern = "evaluation_", replacement = "")
        info <- strsplit(info, split = "_")
        k <- as.integer(info[[1]][2])
        method <- info[[1]][1]
        cat("\nMethod: ", method, "\n\n")
        
        results <- results |> rbind(
          read_tsv(ef, col_names = TRUE) |>
            mutate(VARIANT = factor(VARIANT, levels = variant_regions$VARIANT)) |>
            complete(VARIANT) |>
            mutate(
              COVERAGE = cov,
              TECHNOLOGY = technology,
              SPANNING = spanning,
              ALN_TYPE = aln_type,
              METHOD = method,
              KMER_SIZE = k,
              ACCURACY = acc
            )
        )
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
    # METHOD != "ITER_POA_PERMISSIVE"
  ) %>%
  filter(
    ! (METHOD == "ITER_POA" & TECHNOLOGY == "ont"),
    ! (METHOD == "ITER_POA_PERMISSIVE" & TECHNOLOGY == "pacbio")
  ) %>%
  mutate(
    METHOD = ifelse(TECHNOLOGY == "ont" & METHOD == "ITER_POA_PERMISSIVE", "ITER_POA", METHOD)
  )

# calculate alignment match rate and join TYPE and ZYGOSITY
results <- results %>% mutate(
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

results <- results %>%
  mutate(TECHNOLOGY = toupper(TECHNOLOGY))
results[which(results$TECHNOLOGY == "PACBIO"), "TECHNOLOGY"] <- "PacBio HiFi"

results <- results %>%
  mutate(METHOD = rename_methods(METHOD))

results$ALLELE_MATCH[results$MATCH_RATE < 0.95] <- -1
write_rds(results, paste0(results_dir_kmeans, "/kmer_results.rds"))


grouped_variants <- results %>%
  group_by(
    VARIANT, METHOD, COVERAGE, TECHNOLOGY, HOMOZYGOUS, ALN_TYPE, SPANNING, KMER_SIZE, ACCURACY
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
write_rds(grouped_variants, paste0(results_dir_kmeans, "/kmer_grouped_variants.rds"))

pr_values_het <- grouped_variants %>%
  filter(!HOMOZYGOUS) %>%
  group_by(METHOD, GROUP, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, ACCURACY, KMER_SIZE) %>%
  summarize(N = n(), .groups = "keep") %>%
  pivot_wider(names_from = GROUP, values_from = N, names_prefix = "Group_") %>%
  replace_na(list(Group_A = 0, Group_B = 0, Group_C = 0, Group_D = 0, Group_E = 0, Group_F = 0, Group_G = 0)) %>%
  group_by(METHOD, TECHNOLOGY, ALN_TYPE, SPANNING, COVERAGE, ACCURACY, KMER_SIZE) %>%
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
  mutate(ORIGIN = "kmeans_sim")
pr_values_het %>% write_rds(paste0(results_dir_kmeans, "/pr_values_kmeans.rds"))


true_positives <- results %>%
  filter(!is.na(MATCH_RATE), !is.na(BP_COVERED), BP_COVERED == 1) %>%
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, ALLELE_MATCH, KMER_SIZE, ACCURACY) |>
  summarize(
    ErrorRate = min(1 - MATCH_RATE), 
    .groups = "keep"
  ) |>
  group_by(VARIANT, METHOD, COVERAGE, TYPE, HOMOZYGOUS, TECHNOLOGY, ALN_TYPE, SPANNING, KMER_SIZE, ACCURACY) |>
  summarize(
    ErrorRate = ifelse(HOMOZYGOUS[1] && n() == 2, min(ErrorRate), ErrorRate),
    .groups = "keep"
  )
true_positives %>% write_rds(paste0(results_dir_kmeans, "/true_positives.rds"))
