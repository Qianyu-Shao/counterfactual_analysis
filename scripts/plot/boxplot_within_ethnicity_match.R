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

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
})

# -----------------------------------------
# 0. Configure source pipeline
# -----------------------------------------
source_script <- "pa vs sns.R"
output_name <- if (source_script == "dai_with_temple.R") {
  "within_ethnicity_matched_boxplot_dai_temple.png"
} else {
  "within_ethnicity_matched_boxplot_bulang.png"
}

source("./scripts/function.R")
source(file.path("./scripts", source_script))
source("./scripts/within_ethnicity_match.R")
source("./scripts/extract_encroachment_after_match.R")

# -----------------------------------------
# 1. Build plotting data from matched sample
# -----------------------------------------
make_ethnic_long <- function(forest_df, woody_df, ethnic_name) {
  forest_df <- forest_df %>%
    filter(ethnicity == ethnic_name)

  woody_df <- woody_df %>%
    filter(ethnicity == ethnic_name)

  forest_plot <- forest_df %>%
    st_drop_geometry() %>%
    select(is_sns, loss_1_34_intensity_14, loss_1_345_intensity_14) %>%
    pivot_longer(
      cols = c(loss_1_34_intensity_14, loss_1_345_intensity_14),
      names_to = "outcome",
      values_to = "encroachment"
    )

  woody_plot <- woody_df %>%
    st_drop_geometry() %>%
    select(is_sns, loss_12_34_intensity_14, loss_12_345_intensity_14) %>%
    pivot_longer(
      cols = c(loss_12_34_intensity_14, loss_12_345_intensity_14),
      names_to = "outcome",
      values_to = "encroachment"
    )

  bind_rows(forest_plot, woody_plot) %>%
    mutate(ethnicity = ethnic_name)
}

plot_data <- bind_rows(
  make_ethnic_long(matched_forest_dai_hani, matched_woody_dai_hani, "Dai"),
  make_ethnic_long(matched_forest_dai_hani, matched_woody_dai_hani, "Hani")
) %>%
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
  filter(!is.na(encroachment), !is.na(ethnicity), !is.na(is_sns))

# -----------------------------------------
# 2. P-values and direction arrows
# -----------------------------------------
stat_test <- plot_data %>%
  group_by(ethnicity, outcome) %>%
  wilcox_test(encroachment ~ is_sns) %>%
  ungroup()

direction_df <- plot_data %>%
  group_by(ethnicity, outcome, is_sns) %>%
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
  group_by(ethnicity, outcome) %>%
  summarise(
    y.position = max(encroachment, na.rm = TRUE) * 1.18 + 1e-6,
    .groups = "drop"
  )

stat_test <- stat_test %>%
  left_join(direction_df, by = c("ethnicity", "outcome")) %>%
  left_join(y_pos, by = c("ethnicity", "outcome")) %>%
  mutate(
    xmin = 1,
    xmax = 2,
    p.label = case_when(
      p < 0.001 ~ paste0(direction, " p < 0.001"),
      TRUE ~ paste0(direction, " p = ", formatC(p, format = "f", digits = 3))
    )
  )

# -----------------------------------------
# 3. Plot
# -----------------------------------------
p <- ggplot(plot_data, aes(x = is_sns, y = encroachment, fill = is_sns, color = is_sns)) +
  geom_boxplot(alpha = 0.35, width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.35, size = 1.1) +
  stat_pvalue_manual(
    stat_test,
    label = "p.label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    step.increase = 0,
    inherit.aes = FALSE
  ) +
  facet_grid(ethnicity ~ outcome, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
  scale_color_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
  labs(
    title = "Encroachment Intensity After Within-Ethnicity Matching",
    subtitle = "\u2191 True SNS > Pseudo-SNS; \u2193 True SNS < Pseudo-SNS",
    x = NULL,
    y = "Encroachment Intensity",
    fill = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(
  file.path("./plot", output_name),
  p,
  width = 16,
  height = 10,
  dpi = 300
)
