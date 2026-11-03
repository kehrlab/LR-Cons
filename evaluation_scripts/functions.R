rename_methods <- function(method_names)
{
  method_names[method_names == "KMEANS"] <- "k-mer clust"
  method_names[method_names == "TRUE_HAP"] <- "True phase"
  method_names[method_names == "ITER_POA"] <- "Multi POA"
  method_names[method_names == "fuzzy"] <- "RapidFuzz"
  method_names[method_names == "Molephase"] <- "A priori"
  method_names[method_names == "ITER_POA_PERMISSIVE"] <- "Multi POA v2"
  return(method_names)
}


zygosity_from_gt <- function(gts)
{
  zygosities <- rep("HET", length(gts))
  zygosities[gts == "1|1"] <- "HOM"
  return(zygosities)
}


determine_group <- function(HOMOZYGOUS, N, N_unique, NoSeq, ALL_UNCOVERED, ONE_UNCOVERED, ALL_UNMATCHED, ONE_UNMATCHED)
{
  groups <- rep(NA, length(HOMOZYGOUS))
  
  groups[!is.na(HOMOZYGOUS) & N == 1 & !NoSeq & !HOMOZYGOUS] <- "E"
  groups[!is.na(HOMOZYGOUS) & N == 1 & !NoSeq & HOMOZYGOUS] <- "A"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & HOMOZYGOUS] <- "B"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & N_unique == 1 & !HOMOZYGOUS] <- "F"
  groups[!is.na(HOMOZYGOUS) & N == 2 & !NoSeq & N_unique == 2 & !HOMOZYGOUS] <- "D"
  
  groups[!is.na(HOMOZYGOUS) & N == 1 & NoSeq & !HOMOZYGOUS] <- "G"
  groups[!is.na(HOMOZYGOUS) & N == 1 & NoSeq & HOMOZYGOUS] <- "C"
  groups[!is.na(HOMOZYGOUS) & HOMOZYGOUS & (ALL_UNCOVERED | ALL_UNMATCHED)] <- "C"
  groups[!is.na(HOMOZYGOUS) & !HOMOZYGOUS & (ALL_UNCOVERED | ALL_UNMATCHED)] <- "G"
  
  groups[!is.na(HOMOZYGOUS) &!HOMOZYGOUS & N == 2 & !ALL_UNCOVERED & !ALL_UNMATCHED & (ONE_UNCOVERED | ONE_UNMATCHED)] <- "F"
  
  return(groups)
}