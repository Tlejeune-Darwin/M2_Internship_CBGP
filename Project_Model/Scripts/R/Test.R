library(tidyverse)

# Read the file
df <- read_csv("summary_table.csv")

# Keep interesting columns only
df <- df %>%
  select(simulation_id, census_N, LD_Ne_0.05_Pop1, LD_Ne_0.05_Pop2,
         HE_Neb_mean_Pop1, HE_Neb_mean_Pop2, Coan_Neb_n_Pop1,
         Coan_Neb_n_Pop2, Ne_Pollak, Ne_Nei, Ne_Jorde, mean_alleles_pop1, mean_alleles_pop2)

# True value is important for the bias
true_value <- df$pop_size

# Stat calculation
summary_list <- lapply(df, function(col) {
  if (is.numeric(col)) {
    n <- sum(!is.na(col))
    moy <- mean(col, na.rm = TRUE)
    var_ <- var(col, na.rm = TRUE)
    sd_ <- sd(col, na.rm = TRUE)
    standard_err <- sd_ / sqrt(n)
    t_val <- qt(0.975, df = n - 1)
    ic_lower <- moy - t_val * standard_err
    ic_upper <- moy + t_val * standard_err
    biais <- moy - true_value
    biais_pct <- 100 * biais / true_value
    rmse <- sqrt(mean((col - true_value)^2, na.rm = TRUE))
    rmse_relatif <- rmse / true_value
    min_val <- min(col, na.rm = TRUE)
    max_val <- max(col, na.rm = TRUE)
    
    # Wang formulas
    inv_est <- 1 / col
    inv_est <- inv_est[is.finite(inv_est)]
    inv_mean <- mean(inv_est, na.rm = TRUE)
    Ne_harmonic <- 1 / inv_mean
    bias_wang <- inv_mean - (1 / true_value)
    var_wang <- mean((inv_mean - inv_est)^2, na.rm = TRUE)
    rmse_wang <- sqrt(mean((inv_est - (1 / true_value))^2, na.rm = TRUE))
    
    # Create the summary table
    return(data.frame(
      Moyenne = round(moy, 5),
      Min = round(min_val, 5),
      Max = round(max_val, 5),
      IC_inf = round(ic_lower, 5),
      IC_sup = round(ic_upper, 5),
      Biais = round(biais, 5),
      Biais_pct = paste0(round(biais_pct, 2), "%"),
      RMSE = round(rmse_relatif, 5),
      Wang_Biais = round(bias_wang, 5),
      Wang_Variance = round(var_wang, 5),
      Wang_RMSE = round(rmse_wang, 5)
    ))
  } else {
    return(data.frame(
      Moyenne = NA, Min = NA, Max = NA,
      IC_inf = NA, IC_sup = NA, Biais = NA,
      Biais_pct = NA, RMSE = NA,
      Wang_Biais = NA, Wang_Variance = NA, Wang_RMSE = NA
    ))
  }
})


# Merging of the different rows
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- names(df)

# Print of the table
print(summary_df)


# Sauvegarde du tableau en .csv
# write_csv(summary_df, "Ne_31250.csv")

#######################################################################################

# Chargement des packages
library(tidyverse)

# Lecture du fichier
df <- read_csv("summary_table.csv")

# Calcul du biais
df <- df %>%
  mutate(
    bias = census_N - pop_size,
    relative_bias = (census_N - pop_size) / pop_size
  )

# Statistiques descriptives
summary(df$bias)
summary(df$relative_bias)

# Histogramme du biais absolu
ggplot(df, aes(x = bias)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais absolu (census_N - pop_size)",
       x = "Biais", y = "Nombre de simulations")

# Boxplot du biais relatif
ggplot(df, aes(y = relative_bias)) +
  geom_boxplot(fill = "orange") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Biais relatif de l'estimation du census_N",
       y = "(census_N - pop_size) / pop_size")

# Affichage de quelques exemples
head(df %>% select(simulation_id, pop_size, census_N, bias, relative_bias), 10)

#######################################################################################

# Chargement des packages
library(tidyverse)

# Lecture du fichier
df <- read_csv("summary_table.csv")

# Calcul des biais en excluant les NA
df <- df %>%
  mutate(
    bias_LD_pop1 = LD_Ne_0.05_Pop1 - pop_size,
    relative_bias_LD_pop1 = (LD_Ne_0.05_Pop1 - pop_size) / pop_size,
    
    bias_LD_pop2 = LD_Ne_0.05_Pop2 - pop_size,
    relative_bias_LD_pop2 = (LD_Ne_0.05_Pop2 - pop_size) / pop_size
  )

# Statistiques descriptives avec exclusion des NA
summary(df$bias_LD_pop1)
summary(df$relative_bias_LD_pop1)
summary(df$bias_LD_pop2)
summary(df$relative_bias_LD_pop2)

# Histogrammes avec exclusion des NA
ggplot(df %>% filter(!is.na(bias_LD_pop1)), aes(x = bias_LD_pop1)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Biais absolu LD Ne (pop 1)", x = "Biais", y = "Nombre de simulations")

ggplot(df %>% filter(!is.na(bias_LD_pop2)), aes(x = bias_LD_pop2)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Biais absolu LD Ne (pop 2)", x = "Biais", y = "Nombre de simulations")

# Préparation des biais relatifs en long format (en filtrant les NA)
df_long <- df %>%
  select(simulation_id, pop_size,
         relative_bias_LD_pop1, relative_bias_LD_pop2) %>%
  pivot_longer(cols = starts_with("relative_bias"), names_to = "population", values_to = "relative_bias") %>%
  mutate(population = recode(population,
                             "relative_bias_LD_pop1" = "Pop 1",
                             "relative_bias_LD_pop2" = "Pop 2")) %>%
  filter(!is.na(relative_bias))

# Boxplot
ggplot(df_long, aes(x = population, y = relative_bias)) +
  geom_boxplot(fill = "orange") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Biais relatif LD Ne 0.05",
       y = "(LD_Ne - pop_size) / pop_size", x = "Population")

# Aperçu de quelques lignes complètes
df %>%
  select(simulation_id, pop_size, LD_Ne_0.05_Pop1, LD_Ne_0.05_Pop2,
         bias_LD_pop1, relative_bias_LD_pop1,
         bias_LD_pop2, relative_bias_LD_pop2) %>%
  drop_na() %>%
  head(10)

###################################################################################

library(tidyverse)

# Lecture du fichier CSV
df <- read_csv("summary_table.csv")

# Calcul des biais
df <- df %>%
  mutate(
    bias_pollak = Ne_Pollak - pop_size,
    relative_bias_pollak = (Ne_Pollak - pop_size) / pop_size,
    
    bias_nei = Ne_Nei - pop_size,
    relative_bias_nei = (Ne_Nei - pop_size) / pop_size,
    
    bias_jorde = Ne_Jorde - pop_size,
    relative_bias_jorde = (Ne_Jorde - pop_size) / pop_size
  )

# Pollak
ggplot(df %>% filter(!is.na(bias_pollak)), aes(x = bias_pollak)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais absolu - Méthode Pollak", x = "Biais (Ne_Pollak - pop_size)", y = "Nombre de simulations")

# Nei
ggplot(df %>% filter(!is.na(bias_nei)), aes(x = bias_nei)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais absolu - Méthode Nei", x = "Biais (Ne_Nei - pop_size)", y = "Nombre de simulations")

# Jorde-Ryman
ggplot(df %>% filter(!is.na(bias_jorde)), aes(x = bias_jorde)) +
  geom_histogram(bins = 30, fill = "orange", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais absolu - Méthode Jorde-Ryman", x = "Biais (Ne_Jorde - pop_size)", y = "Nombre de simulations")

summary(select(df, bias_pollak, relative_bias_pollak,
               bias_nei, relative_bias_nei,
               bias_jorde, relative_bias_jorde))

# Histogramme du biais relatif - Pollak
ggplot(df %>% filter(!is.na(relative_bias_pollak)), aes(x = relative_bias_pollak)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais relatif - Méthode Pollak",
       x = "(Ne_Pollak - pop_size) / pop_size", y = "Nombre de simulations")

# Histogramme du biais relatif - Nei
ggplot(df %>% filter(!is.na(relative_bias_nei)), aes(x = relative_bias_nei)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais relatif - Méthode Nei",
       x = "(Ne_Nei - pop_size) / pop_size", y = "Nombre de simulations")

# Histogramme du biais relatif - Jorde-Ryman
ggplot(df %>% filter(!is.na(relative_bias_jorde)), aes(x = relative_bias_jorde)) +
  geom_histogram(bins = 30, fill = "orange", color = "black") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Biais relatif - Méthode Jorde-Ryman",
       x = "(Ne_Jorde - pop_size) / pop_size", y = "Nombre de simulations")

df_selection <- df %>%
  select(simulation_id, pop_size, census_N, Ne_Jorde, Ne_Nei, Ne_Pollak, relative_bias_jorde, relative_bias_nei, relative_bias_pollak)

ggplot(df_long, aes(x = relative_bias, fill = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution du biais relatif (densité)",
       x = "Biais relatif", fill = "Méthode")

########################################################################################

library(tidyverse)

# Lecture du fichier
df <- read_csv("summary_table.csv")

# Calcul des biais relatifs si besoin
df <- df %>%
  mutate(
    relative_bias_pollak = (Ne_Pollak - pop_size) / pop_size,
    relative_bias_nei    = (Ne_Nei    - pop_size) / pop_size,
    relative_bias_jorde  = (Ne_Jorde  - pop_size) / pop_size
  )

# Mise en format long + suppression des NA
df_long <- df %>%
  select(simulation_id, pop_size,
         relative_bias_pollak, relative_bias_nei, relative_bias_jorde) %>%
  pivot_longer(cols = starts_with("relative_bias"),
               names_to = "method", values_to = "relative_bias") %>%
  mutate(method = recode(method,
                         "relative_bias_pollak" = "Pollak",
                         "relative_bias_nei" = "Nei",
                         "relative_bias_jorde" = "Jorde")) %>%
  drop_na(relative_bias)
df_long <- df_long %>%
  filter(simulation_id != "sim_20250416_091816_local")%>%
  filter(simulation_id != "sim_20250416_091919_local")


# Scatter plot
ggplot(df_long, aes(x = pop_size, y = relative_bias, color = method)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10() +
  labs(title = "Biais relatif des méthodes temporelles vs pop_size",
       x = "Taille de population simulée (log scale)",
       y = "Biais relatif", color = "Méthode")
