install.packages("abcrf", dependencies = TRUE)
library(tidyverse)
library(abcrf)

# Read the file
ref_data <- read_csv("Bureau/simulations/Ref_table/summary_table.csv")

params <- ref_data[, c("simulation_id", "pop_size", "num_loci", "sample1_size_Ne", "sample2_size_Ne", "sample1_size_CMR", "sample2_size_CMR", "mutation_rate", "recap_Ne")]
stats_keywords <- c("id", "LD", "HE", "Coan", "het", "alleles", "P", "N", "J")
is_name_col <- sapply(names(ref_data), function(col) {
  any(sapply(stats_keywords, function(kw) grepl(kw, col)))
})
stats_table <- ref_data[, is_name_col]
stats_table$census_N[stats_table$census_N == 0] <- NA
