install.packages("abcrf", dependencies = TRUE)
library(tidyverse)
library(abcrf)

# Read the file
ref_data <- read_csv("Bureau/simulations/Ref_table/summary_table.csv")

params <- ref_data[, c("simulation_id", "pop_size", "num_loci", "sample1_size_Ne", "sample2_size_Ne", "sample1_size_CMR", "sample2_size_CMR", "mutation_rate", "recap_Ne")]
stats_keywords <- c("id", "LD", "HE", "Coan", "het", "alleles", "P", "N", "J")
stats_keywords <- c("id", "het", "alleles", "P", "N", "J")
is_name_col <- sapply(names(ref_data), function(col) {
  any(sapply(stats_keywords, function(kw) grepl(kw, col)))
})
stats_table <- ref_data[, is_name_col]
stats_table$census_N[stats_table$census_N == 0] <- NA
reference_imputed <- stats_table  # copie de travail

for(col in names(reference_imputed)) {
  if(anyNA(reference_imputed[[col]])) {
    med <- median(reference_imputed[[col]], na.rm = TRUE)
    reference_imputed[[col]][is.na(reference_imputed[[col]])] <- med
  }
}

# ABCRF

target_param <- "pop_size"
learning_data <- data.frame(y = params[[target_param]], stats_table)
names(learning_data)[1] <- target_param

formula_rf <- as.formula(paste(target_param, "~ ."))
model_rf <- regAbcrf(formula = formula_rf, data = learning_data, ntree = 500)
summary(model_rf)
length(model.response(mf))                 # Vecteur des paramètres simulés
length(model.rf$predictions)               # Vecteur des prédictions (OOB)
# Supposons que ton tableau s'appelle "reference"
na_counts <- colSums(is.na(learning_data))

# Afficher uniquement les colonnes avec au moins un NA
na_columns <- na_counts[na_counts > 0]

# Résultat lisible
if(length(na_columns) == 0) {
  cat("✅ Aucune colonne ne contient de NA.\n")
} else {
  cat("⚠️ Colonnes contenant des NA :\n")
  print(na_columns)
}
