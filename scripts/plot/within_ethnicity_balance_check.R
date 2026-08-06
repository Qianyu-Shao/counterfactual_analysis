library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readr)

calc_smd <- function(data, var, treat = "is_sns") {
  x1 <- data[[var]][data[[treat]] == 1]
  x0 <- data[[var]][data[[treat]] == 0]
  
  m1 <- mean(x1, na.rm = TRUE)
  m0 <- mean(x0, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE)
  s0 <- sd(x0, na.rm = TRUE)
  
  pooled_sd <- sqrt((s1^2 + s0^2) / 2)
  smd <- ifelse(is.na(pooled_sd) || pooled_sd == 0, NA, (m1 - m0) / pooled_sd)
  
  tibble(
    variable = var,
    mean_true = m1,
    mean_pseudo = m0,
    smd = smd
  )
}

make_ethnicity <- function(df) {
  df %>%
    mutate(
      ethnicity = case_when(
        dai == 1 ~ "Dai",
        hani == 1 ~ "Hani",
        TRUE ~ "Bulang"
      )
    )
}

check_balance_by_ethnicity <- function(df, vars, sample_name) {
  df %>%
    group_split(ethnicity) %>%
    map_dfr(function(subdf) {
      eth <- unique(subdf$ethnicity)
      map_dfr(vars, ~calc_smd(subdf, .x)) %>%
        mutate(
          ethnicity = eth,
          sample = sample_name
        )
    }) %>%
    select(sample, ethnicity, everything())
}

plot_balance_density <- function(df, vars, sample_name, file_name) {
  plot_data <- df %>%
    select(all_of(vars), is_sns, ethnicity) %>%
    pivot_longer(
      cols = all_of(vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(!is.na(value), !is.na(is_sns)) %>%
    mutate(
      is_sns = factor(is_sns, levels = c(0, 1), labels = c("Pseudo-SNS", "True SNS")),
      variable = case_when(
        variable == "mean_elev" ~ "Elevation",
        variable == "mean_slope" ~ "Slope",
        variable == "dist_to_road" ~ "Distance to Road",
        variable == "mode_aspect" ~ "Aspect",
        variable == "forest_area_02" ~ "Initial Forest Area",
        variable == "woody_area_02" ~ "Initial Woody Area",
        variable == "dist_to_village" ~ "Distance to Village",
        TRUE ~ variable
      ),
      value_plot = case_when(
        variable %in% c("Distance to Road", "Distance to Village", "Initial Forest Area", "Initial Woody Area") ~ log1p(value),
        TRUE ~ value
      )
    )
  
  p <- ggplot(plot_data, aes(x = value_plot, color = is_sns, fill = is_sns)) +
    geom_density(alpha = 0.2, adjust = 1.1) +
    facet_grid(ethnicity ~ variable, scales = "free") +
    scale_color_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
    scale_fill_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
    labs(
      title = paste("Within-Ethnicity Covariate Balance:", sample_name),
      subtitle = "Distance and area variables shown on log1p scale where needed",
      x = NULL,
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(file_name, p, width = 14, height = 9, dpi = 300)
  
  p
}

# ----------------------------
# 1. Prepare matched datasets
# ----------------------------

forest_df <- matched_data_12345_forest_500 %>%
  sf::st_drop_geometry() %>%
  make_ethnicity()

woody_df <- matched_data_12345_woody_500 %>%
  sf::st_drop_geometry() %>%
  make_ethnicity()

forest_vars <- c(
  "mean_elev",
  "mean_slope",
  "dist_to_road",
  "mode_aspect",
  "forest_area_02",
  "dist_to_village"
)

woody_vars <- c(
  "mean_elev",
  "mean_slope",
  "dist_to_road",
  "mode_aspect",
  "woody_area_02",
  "dist_to_village"
)

# ----------------------------
# 2. Balance tables
# ----------------------------

forest_balance <- check_balance_by_ethnicity(forest_df, forest_vars, "forest matched sample")
woody_balance  <- check_balance_by_ethnicity(woody_df, woody_vars, "woody matched sample")

balance_table <- bind_rows(forest_balance, woody_balance) %>%
  mutate(
    abs_smd = abs(smd),
    balance_flag = case_when(
      abs_smd < 0.1 ~ "Good",
      abs_smd < 0.2 ~ "Acceptable",
      TRUE ~ "Imbalanced"
    )
  )

print(balance_table)

# Save tables
write_csv(balance_table, "./plot/within_ethnicity_balance_table.csv")
#write_csv(forest_balance, "./plot/within_ethnicity_balance_forest.csv")
#write_csv(woody_balance, "./plot/within_ethnicity_balance_woody.csv")

# ----------------------------
# 3. Density plots
# ----------------------------

plot_balance_density(
  forest_df,
  forest_vars,
  "Forest matched sample",
  "./plot/within_ethnicity_balance_forest.png"
)

plot_balance_density(
  woody_df,
  woody_vars,
  "Woody matched sample",
  "./plot/within_ethnicity_balance_woody.png"
)
