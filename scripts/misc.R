library(tidyverse)

my_colors <- c(
  '#ffd266ff', # gelb
  '#f58080ff', #rot
  '#6a00e5ff',
  '#12beccff',
  '#fcab79ff',
  '#0059ecff',
  '#000c52',
  "lightgreen", # not LIT
  '#ff69b4' # not LIT
)

method_order <- c(
  "abPOA",
  "Multi POA",
  "RapidFuzz",
  "k-mer clust",
  "A priori",
  "True phase"
)

method_order_old <- c(
  "TRUE_HAP",
  "ITER_POA",
  "fuzzy",
  "Molephase",
  "KMEANS",
  "abPOA"
)



group_order <- c(4, 5, 2, 7, 3, 6, 1)
group_names <- c(
  "HOM_ONE", "HOM_TWO", "HOM_NONE",
  "HET_DISTINCT", "HET_ONE", "HET_REDUNDANT", "HET_NONE"
)

custom_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "white")
  )

rotate_text_x <- theme(
  axis.text.x = element_text(angle = 20, hjust = 0.9, vjust = 0.9),
  axis.title.x = element_blank()
)

empty_x <- theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank()
)
