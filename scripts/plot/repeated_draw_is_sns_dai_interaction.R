library(dplyr)
library(ggplot2)
library(readr)

base_plot_dir <- "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis/plot"
out_dir <- file.path(base_plot_dir, "repeated_draw_is_sns_dai_interaction")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source_dirs <- list(
  bulang = file.path(base_plot_dir, "robust_check_50_m3_slm_bulang"),
  altcoding = file.path(base_plot_dir, "robust_check_50_m3_slm_bulang_altcoding")
)

coef_data <- bind_rows(lapply(names(source_dirs), function(src_name) {
  read_csv(file.path(source_dirs[[src_name]], "draw_level_coefficients.csv"), show_col_types = FALSE) %>%
    mutate(source_version = src_name)
}))

spec_df <- tibble(
  source_version = c("altcoding", "altcoding", "bulang", "bulang"),
  sample = c("Forest", "Forest", "Woody", "Woody"),
  outcome = c(
    "forest_to_cashcrop",
    "forest_to_farmland",
    "woody_to_cashcrop",
    "woody_to_farmland"
  ),
  term = "is_sns:dai_eth",
  x = 2:5,
  label = c(
    "Forest to\nrubber",
    "Forest to\nfarmland",
    "Shrubland to\nrubber",
    "Shrubland to\nfarmland"
  )
)

effect_df <- bind_rows(lapply(seq_len(nrow(spec_df)), function(i) {
  spec_row <- spec_df[i, ]
  draw_df <- coef_data %>%
    filter(
      source_version == spec_row$source_version,
      sample == spec_row$sample,
      outcome == spec_row$outcome,
      term == spec_row$term
    )

  tibble(
    x = spec_row$x,
    label = spec_row$label,
    mean_estimate = mean(draw_df$estimate, na.rm = TRUE) * 100,
    ci_low = quantile(draw_df$estimate, 0.025, na.rm = TRUE) * 100,
    ci_high = quantile(draw_df$estimate, 0.975, na.rm = TRUE) * 100,
    p_ratio = mean(draw_df$p_value < 0.05, na.rm = TRUE),
    p_label = sprintf("%.2f", mean(draw_df$p_value < 0.05, na.rm = TRUE))
  )
}))

baseline_df <- tibble(
  x = 1,
  label = "Hani",
  mean_estimate = 0,
  ci_low = 0,
  ci_high = 0,
  p_ratio = NA_real_,
  p_label = NA_character_
)

plot_df <- bind_rows(baseline_df, effect_df)

y_min <- min(plot_df$ci_low, na.rm = TRUE)
y_max <- max(plot_df$ci_high, na.rm = TRUE)
y_pad <- max((y_max - y_min) * 0.18, 4)
y_top <- y_max + y_pad * 1.15

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.8) +
  geom_linerange(
    data = effect_df,
    aes(x = x, ymin = ci_low, ymax = ci_high),
    linewidth = 0.8,
    color = "grey45"
  ) +
  geom_point(
    data = baseline_df,
    aes(x = x, y = mean_estimate),
    shape = 21,
    size = 3.6,
    stroke = 0.9,
    fill = "white",
    color = "grey25"
  ) +
  geom_point(
    data = effect_df,
    aes(x = x, y = mean_estimate),
    shape = 19,
    size = 4.2,
    color = "grey20"
  ) +
  geom_text(
    data = effect_df,
    aes(x = x, y = y_top, label = p_label),
    size = 5.2,
    fontface = "bold",
    vjust = 0
  ) +
  scale_x_continuous(
    breaks = plot_df$x,
    labels = plot_df$label,
    limits = c(0.85, 5.15),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(y_min - y_pad * 0.25, y_top + y_pad * 0.95)
  ) +
  labs(
    x = NULL,
    y = "Encroachment Intensity Difference (%)",
    title = "Dai-Hani Sacred-site Effect Difference Across Land-conversion Types"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.line.x = element_line(color = "grey35", linewidth = 0.7),
    axis.ticks.x = element_line(color = "grey35", linewidth = 0.7),
    axis.ticks.length.x = unit(5, "pt"),
    plot.margin = margin(12, 16, 10, 12)
  )

ggsave(
  filename = file.path(out_dir, "is_sns_dai_interaction_summary.png"),
  plot = p,
  width = 10.4,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "is_sns_dai_interaction_summary.pdf"),
  plot = p,
  width = 10.4,
  height = 5.2,
  bg = "white"
)

write_csv(plot_df, file.path(out_dir, "is_sns_dai_interaction_summary_data.csv"))

message("Saved outputs to: ", out_dir)
