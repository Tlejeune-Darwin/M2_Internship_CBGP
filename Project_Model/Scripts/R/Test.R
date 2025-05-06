install.packages("abcrf", dependencies = TRUE)
library(tidyverse)
library(abcrf)

# Read the file
ref_data <- read_csv("Bureau/simulations/summary_table_all_batches.csv")

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

#########################################

data <- read.csv("Bureau/simulations/Ref_table/summary_table.csv", header = TRUE)

# 2. Définir les colonnes à supprimer
colonnes_a_supprimer <- c(
  "HE_Neb_mean_0.050_Pop2",
  "HE_Neb_mean_0.020_Pop2",
  "HE_Neb_mean_0.010_Pop2",
  "HE_Neb_mean_0.000_Pop2",
  "HE_Neb_mean_0.050_Pop1",
  "HE_Neb_mean_0.020_Pop1",
  "HE_Neb_mean_0.010_Pop1",
  "HE_Neb_mean_0.000_Pop1",
  "Coan_Neb_n_Pop1",
  "Coan_Neb_n_Pop2",
  "LD_Ne_0.050_Pop1",
  "LD_Ne_0.020_Pop1",
  "LD_Ne_0.010_Pop1",
  "LD_Ne_0.000_Pop1",
  "LD_Ne_0.050_Pop2",
  "LD_Ne_0.020_Pop2",
  "LD_Ne_0.010_Pop2",
  "LD_Ne_0.000_Pop2"
)

# 3. Supprimer ces colonnes
data_clean <- data[, !(names(data) %in% colonnes_a_supprimer)]
data_clean <- na.omit(data_clean)

params <- data_clean[, c("simulation_id", "pop_size", "num_loci", "sample1_size_Ne", "sample2_size_Ne", "sample1_size_CMR", "sample2_size_CMR", "mutation_rate", "recap_Ne")]
stats_keywords <- c("id", "LD", "HE", "Coan", "het", "alleles", "P_", "N_", "J_")
#stats_keywords <- c("id", "het", "alleles", "P", "N", "J")
is_name_col <- sapply(names(data_clean), function(col) {
  any(sapply(stats_keywords, function(kw) grepl(kw, col)))
})
stats_table <- data_clean[, is_name_col]
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
model_rf$model.rf$prediction.error
importance_stats <- sort(model_rf$model.rf$variable.importance, decreasing = T)
importance_df <- data.frame(
  Statistic = names(importance_stats),
  Importance = as.numeric(importance_stats)
)
head(importance_stats, 30)
mean(learning_data$pop_size)
var(learning_data$pop_size)

# 1. Prédire en utilisant le modèle sur les mêmes données
predictions <- predict(model_rf, training = learning_data, obs = learning_data)$expectation
head(predictions)
plot(learning_data$pop_size, predictions,
     xlab = "Vraie pop_size", ylab = "Prédiction pop_size (OOB)",
     main = "Comparaison vraie vs prédite",
     pch = 20, col = "blue")
abline(0, 1, col = "red")

# Afficher uniquement les colonnes avec au moins un NA
na_columns <- na_counts[na_counts > 0]

# Résultat lisible
if(length(na_columns) == 0) {
  cat("✅ Aucune colonne ne contient de NA.\n")
} else {
  cat("⚠️ Colonnes contenant des NA :\n")
  print(na_columns)
}

