`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
}

if (!is.null(script_file) && nzchar(script_file)) {
  script_file <- gsub("~\\+~", " ", script_file, fixed = FALSE)
}

project_root <- if (!is.null(script_file) && nzchar(script_file)) {
  normalizePath(
    file.path(dirname(script_file), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  normalizePath(
    "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis",
    winslash = "/",
    mustWork = TRUE
  )
}

setwd(project_root)

invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
})

# -----------------------------------------
# 0. Configure which pseudo-site pipeline
# -----------------------------------------
# Change this to "bulang_clean.R" if needed.
source_script <- "bulang_clean.R"
plot_seed <- 20260510
output_name <- if (source_script == "dai_with_temple.R") {
  "overlap_ethnicity_boxplot_p_arrow_dai_temple.png"
} else {
  "overlap_ethnicity_boxplot_p_arrow_bulang.png"
}

set.seed(plot_seed)
message("Using fixed pseudo-site seed: ", plot_seed)

source("./scripts/function.R")
source(file.path("./scripts", source_script))
source("./scripts/within_ethnicity_match.R")
source("./scripts/extract_encroachment_after_match.R")

# -----------------------------------------
# 1. Rebuild overlap samples
# -----------------------------------------
prepare_overlap_sample <- function(data) {
  data %>%
    filter(dai == 1 | hani == 1) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

build_ethnicity_overlap <- function(df) {
  base_df <- df %>%
    st_drop_geometry() %>%
    mutate(
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )

  ps_fit <- glm(
    dai_eth ~ mean_elev + mean_slope + aspect_cos + aspect_sin +
      dist_to_road + dist_to_village,
    family = binomial(),
    data = base_df
  )

  base_df$ps_eth <- predict(ps_fit, type = "response")

  dai_range <- range(base_df$ps_eth[base_df$dai_eth == 1], na.rm = TRUE)
  hani_range <- range(base_df$ps_eth[base_df$dai_eth == 0], na.rm = TRUE)

  lower <- max(dai_range[1], hani_range[1])
  upper <- min(dai_range[2], hani_range[2])

  overlap_ids <- base_df %>%
    filter(ps_eth >= lower, ps_eth <= upper) %>%
    pull(ID)

  df %>%
    filter(ID %in% overlap_ids) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

forest_overlap_df <- matched_forest_dai_hani %>%
  prepare_overlap_sample() %>%
  build_ethnicity_overlap()

woody_overlap_df <- matched_woody_dai_hani %>%
  prepare_overlap_sample() %>%
  build_ethnicity_overlap()

# -----------------------------------------
# 1b. Balance check after Dai-Hani overlap
# -----------------------------------------
calc_dai_hani_balance <- function(data, sample_name, stage_name, vars) {
  clean_df <- data %>%
    prepare_overlap_sample() %>%
    st_drop_geometry() %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )

  n_dai <- sum(clean_df$dai_eth == 1, na.rm = TRUE)
  n_hani <- sum(clean_df$hani_eth == 1, na.rm = TRUE)

  bind_rows(lapply(vars, function(var) {
    dai_values <- clean_df[[var]][clean_df$dai_eth == 1]
    hani_values <- clean_df[[var]][clean_df$hani_eth == 1]

    mean_dai <- mean(dai_values, na.rm = TRUE)
    mean_hani <- mean(hani_values, na.rm = TRUE)
    sd_dai <- sd(dai_values, na.rm = TRUE)
    sd_hani <- sd(hani_values, na.rm = TRUE)
    pooled_sd <- sqrt((sd_dai^2 + sd_hani^2) / 2)
    smd <- ifelse(is.na(pooled_sd) || pooled_sd == 0, NA_real_, (mean_dai - mean_hani) / pooled_sd)

    tibble(
      sample = sample_name,
      stage = stage_name,
      covariate = var,
      n_dai = n_dai,
      n_hani = n_hani,
      mean_dai = mean_dai,
      mean_hani = mean_hani,
      smd = smd,
      abs_smd = abs(smd),
      balance_flag = case_when(
        abs_smd < 0.1 ~ "Good",
        abs_smd < 0.2 ~ "Acceptable",
        TRUE ~ "Imbalanced"
      )
    )
  }))
}

overlap_balance_vars <- c(
  "mean_elev",
  "mean_slope",
  "aspect_cos",
  "aspect_sin",
  "dist_to_road",
  "dist_to_village"
)

overlap_balance_table <- bind_rows(
  calc_dai_hani_balance(matched_forest_dai_hani, "Forest", "Before overlap restriction", overlap_balance_vars),
  calc_dai_hani_balance(forest_overlap_df, "Forest", "After overlap restriction", overlap_balance_vars),
  calc_dai_hani_balance(matched_woody_dai_hani, "Shrubland", "Before overlap restriction", overlap_balance_vars),
  calc_dai_hani_balance(woody_overlap_df, "Shrubland", "After overlap restriction", overlap_balance_vars)
) %>%
  mutate(
    across(c(mean_dai, mean_hani, smd, abs_smd), ~ round(.x, 3))
  )

write.csv(
  overlap_balance_table,
  "./plot/overlap_ethnicity_dai_hani_balance_check.csv",
  row.names = FALSE
)

print(overlap_balance_table)
message("Saved Dai-Hani overlap balance check:")
message("./plot/overlap_ethnicity_dai_hani_balance_check.csv")

# -----------------------------------------
# 2. Build plotting data
# -----------------------------------------
forest_plot <- forest_overlap_df %>%
  st_drop_geometry() %>%
  select(ethnicity, is_sns, loss_1_34_intensity_14, loss_1_345_intensity_14) %>%
  pivot_longer(
    cols = c(loss_1_34_intensity_14, loss_1_345_intensity_14),
    names_to = "outcome",
    values_to = "encroachment"
  )

woody_plot <- woody_overlap_df %>%
  st_drop_geometry() %>%
  select(ethnicity, is_sns, loss_12_34_intensity_14, loss_12_345_intensity_14) %>%
  pivot_longer(
    cols = c(loss_12_34_intensity_14, loss_12_345_intensity_14),
    names_to = "outcome",
    values_to = "encroachment"
  )

plot_data <- bind_rows(forest_plot, woody_plot) %>%
  mutate(
    ethnicity = factor(ethnicity, levels = c("Dai", "Hani")),
    is_sns = factor(is_sns, levels = c(0, 1), labels = c("Pseudo-SNS", "True SNS")),
    outcome = factor(
      outcome,
      levels = c(
        "loss_1_34_intensity_14",
        "loss_1_345_intensity_14",
        "loss_12_34_intensity_14",
        "loss_12_345_intensity_14"
      ),
      labels = c(
        "Forest to Rubber",
        "Forest to Farmland",
        "Shrubland to Rubber",
        "Shrubland to Farmland"
      )
    )
  ) %>%
  filter(!is.na(encroachment), !is.na(ethnicity), !is.na(is_sns)) %>%
  mutate(
    outcome = factor(
      as.character(outcome),
      levels = c(
        "Forest to Rubber",
        "Forest to Farmland",
        "Shrubland to Rubber",
        "Shrubland to Farmland"
      )
    ),
    outcome_panel = factor(
      case_when(
        outcome == "Forest to Rubber" ~ "01 Forest to Rubber",
        outcome == "Forest to Farmland" ~ "02 Forest to Farmland",
        outcome == "Shrubland to Rubber" ~ "03 Shrubland to Rubber",
        outcome == "Shrubland to Farmland" ~ "04 Shrubland to Farmland",
        TRUE ~ NA_character_
      ),
      levels = c(
        "01 Forest to Rubber",
        "02 Forest to Farmland",
        "03 Shrubland to Rubber",
        "04 Shrubland to Farmland"
      )
    ),
    sns_x = ifelse(is_sns == "Pseudo-SNS", 0.88, 1.12)
  )

# -----------------------------------------
# 3. P-values and direction arrows
# -----------------------------------------
stat_test <- plot_data %>%
  group_by(ethnicity, outcome, outcome_panel) %>%
  wilcox_test(encroachment ~ is_sns) %>%
  ungroup()

direction_df <- plot_data %>%
  group_by(ethnicity, outcome, outcome_panel, is_sns) %>%
  summarise(mean_enc = mean(encroachment, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = is_sns, values_from = mean_enc) %>%
  mutate(
    diff = `True SNS` - `Pseudo-SNS`,
    direction = case_when(
      diff > 0 ~ "\u2191",
      diff < 0 ~ "\u2193",
      TRUE ~ "="
    )
  )

y_pos <- plot_data %>%
  group_by(ethnicity, outcome, outcome_panel) %>%
  summarise(
    y_min = min(encroachment, na.rm = TRUE),
    y_max = max(encroachment, na.rm = TRUE),
    y_range = y_max - y_min,
    y.position = y_max + pmax(y_range, 0.05) * 0.12,
    y_n = y_min - pmax(y_range, 0.05) * 0.14,
    .groups = "drop"
  )

stat_test <- stat_test %>%
  left_join(direction_df, by = c("ethnicity", "outcome", "outcome_panel")) %>%
  left_join(y_pos, by = c("ethnicity", "outcome", "outcome_panel")) %>%
  mutate(
    x = 1.00,
    p.label = case_when(
      p < 0.001 ~ paste0(direction, " p < 0.001"),
      TRUE ~ paste0(direction, " p = ", formatC(p, format = "f", digits = 3))
    )
  )

n_labels <- plot_data %>%
  count(ethnicity, outcome, outcome_panel, is_sns, name = "n") %>%
  left_join(y_pos, by = c("ethnicity", "outcome", "outcome_panel")) %>%
  mutate(
    x = ifelse(is_sns == "Pseudo-SNS", 0.88, 1.12),
    label = paste0("n = ", n)
  )

# -----------------------------------------
# 4. Plot
# -----------------------------------------
p <- ggplot(plot_data, aes(x = sns_x, y = encroachment, fill = is_sns)) +
  geom_boxplot(
    width = 0.11,
    linewidth = 0.45,
    color = "grey20",
    outlier.shape = NA,
    alpha = 0.78
  ) +
  geom_jitter(
    aes(color = is_sns),
    width = 0.035,
    alpha = 0.38,
    size = 1.15,
    show.legend = FALSE
  ) +
  geom_text(
    data = stat_test,
    aes(x = x, y = y.position, label = p.label),
    inherit.aes = FALSE,
    size = 5.6,
    fontface = "plain"
  ) +
  geom_text(
    data = n_labels,
    aes(x = x, y = y_n, label = label),
    inherit.aes = FALSE,
    size = 5.0
  ) +
  facet_grid(
    outcome_panel ~ ethnicity,
    scales = "free_y",
    labeller = labeller(outcome_panel = function(x) sub("^\\d+\\s+", "", x))
  ) +
  scale_fill_manual(values = c("Pseudo-SNS" = "#D4D4CC", "True SNS" = "#86AE8E")) +
  scale_color_manual(values = c("Pseudo-SNS" = "#85857E", "True SNS" = "#3F7A4D")) +
  scale_x_continuous(
    breaks = c(0.88, 1.12),
    labels = c("Pseudo-SNS", "True SNS"),
    limits = c(0.72, 1.28)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.20, 0.16))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Encroachment Intensity After Dai-Hani Overlap Restriction",
    subtitle = "\u2191 True SNS > Pseudo-SNS; \u2193 True SNS < Pseudo-SNS",
    x = NULL,
    y = "Encroachment Intensity",
    fill = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.75, "lines"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    plot.margin = margin(12, 18, 22, 12)
  )

ggsave(
  file.path("./plot", output_name),
  p,
  width = 12,
  height = 16,
  dpi = 300
)

zero_share_table <- plot_data %>%
  group_by(ethnicity, outcome, is_sns) %>%
  summarise(
    n = n(),
    n_zero = sum(encroachment == 0, na.rm = TRUE),
    pct_zero = 100 * n_zero / n,
    .groups = "drop"
  ) %>%
  mutate(pct_zero = round(pct_zero, 1))

print(zero_share_table)
