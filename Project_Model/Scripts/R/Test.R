library(tidyverse)

# Read the file
df <- read_csv("simulations/summary_table.csv")

# Keep interesting columns only
df <- df %>%
  select(simulation_id, census_N, LD_Ne_0.05_Pop1, LD_Ne_0.05_Pop2,
         HE_Neb_mean_Pop1, HE_Neb_mean_Pop2, Coan_Neb_n_Pop1,
         Coan_Neb_n_Pop2, Ne_Pollak, Ne_Nei, Ne_Jorde, mean_alleles_pop1, mean_alleles_pop2)

# True value is important for the bias
true_value <- 100

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

library(tidyverse)
library(tibble)

# Chargement manuel de chaque fichier
df_5 <- read_csv("L_5.csv")
df_10 <- read_csv("L_10.csv")
df_15 <- read_csv("L_15.csv")
df_20 <- read_csv("L_20.csv")
df_25 <- read_csv("L_25.csv")
df_30 <- read_csv("L_30.csv")
df_35 <- read_csv("L_35.csv")
df_40 <- read_csv("L_40.csv")

val_5 <- df_5[4, "Wang_RMSE"]
val_10 <- df_10[4, "Wang_RMSE"]
val_15 <- df_15[4, "Wang_RMSE"]
val_20 <- df_20[4, "Wang_RMSE"]
val_25 <- df_25[4, "Wang_RMSE"]
val_30 <- df_30[4, "Wang_RMSE"]
val_35 <- df_35[4, "Wang_RMSE"]
val_40 <- df_40[4, "Wang_RMSE"]


# Création du tableau des résultats
rmse_df <- tibble(
  Param = c(5, 10, 15, 20, 25, 30, 35, 40),
  Wang_RMSE = c(val_5, val_10, val_15, val_20, val_25, val_30, val_35, val_40)
)

library(ggplot2)

# Si ce n’est pas encore fait :
rmse_df$Param <- as.numeric(rmse_df$Param)  # assure que l'axe x est bien numérique
rmse_df$Wang_RMSE <- as.numeric(rmse_df$Wang_RMSE)

library(ggplot2)

ggplot(rmse_df, aes(x = Param, y = Wang_RMSE)) +
  geom_line(linewidth = 1, color = "black", linetype = "dashed") +
  geom_point(size = 2, color = "black") +
  scale_x_log10(breaks = rmse_df$Param) +  # Échelle log avec les valeurs entières affichées
  labs(
    x = expression(italic(L)),
    y = "RMSE"
  ) +
  theme_minimal(base_size = 14)

#######################################################################################

library(tidyverse)
library(tibble)

# Chargement manuel de chaque fichier
df_10 <- read_csv("Ne_10.csv")
df_50 <- read_csv("Ne_50.csv")
df_250 <- read_csv("Ne_250.csv")
df_1250 <- read_csv("Ne_1250.csv")
df_6250 <- read_csv("Ne_6250.csv")
df_31250 <- read_csv("Ne_31250.csv")

val_10 <- df_10[4, "Wang_RMSE"]
val_50 <- df_50[4, "Wang_RMSE"]
val_250 <- df_250[4, "Wang_RMSE"]
val_1250 <- df_1250[4, "Wang_RMSE"]
val_6250 <- df_6250[4, "Wang_RMSE"]
val_31250 <- df_31250[4, "Wang_RMSE"]

# Création du tableau des résultats
rmse_df <- tibble(
  Param = c(10, 50, 250, 1250, 6250, 31250),
  Wang_RMSE = c(val_10, val_50, val_250, val_1250, val_6250, val_31250)
)

library(ggplot2)

# Si ce n’est pas encore fait :
rmse_df$Param <- as.numeric(rmse_df$Param)  # assure que l'axe x est bien numérique
rmse_df$Wang_RMSE <- as.numeric(rmse_df$Wang_RMSE)

library(ggplot2)

ggplot(rmse_df, aes(x = Param, y = Wang_RMSE)) +
  geom_line(linewidth = 1, color = "black", linetype = "dashed") +
  geom_point(size = 2, color = "black") +
  scale_x_log10(breaks = rmse_df$Param) +  # Échelle log avec les valeurs entières affichées
  labs(
    x = expression(italic(Ne)),
    y = "RMSE"
  ) +
  theme_minimal(base_size = 14)

