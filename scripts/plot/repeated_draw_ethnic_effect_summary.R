library(dplyr)
library(ggplot2)
library(readr)

base_plot_dir <- "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis/plot"
out_dir <- file.path(base_plot_dir, "repeated_draw_ethnic_effect_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source_dirs <- list(
  bulang = file.path(base_plot_dir, "robust_check_50_m3_slm_bulang"),
  altcoding = file.path(base_plot_dir, "robust_check_50_m3_slm_bulang_altcoding")
)

coef_data <- bind_rows(lapply(names(source_dirs), function(src_name) {
  read_csv(file.path(source_dirs[[src_name]], "draw_level_coefficients.csv"), show_col_types = FALSE) %>%
    mutate(source_version = src_name)
}))

count_data <- bind_rows(lapply(names(source_dirs), function(src_name) {
  read_csv(file.path(source_dirs[[src_name]], "sample_count_summary.csv"), show_col_types = FALSE) %>%
    mutate(source_version = src_name)
}))

panel_map <- tibble(
  panel_id = c("a", "b", "c", "d"),
  panel_title = c("Dai forest", "Dai shrubland", "Hani forest", "Hani shrubland"),
  source_version = c("bulang", "altcoding", "altcoding", "bulang"),
  sample = c("Forest", "Woody", "Forest", "Woody"),
  true_n_col = c("mean_n_dai_true", "mean_n_dai_true", "mean_n_hani_true", "mean_n_hani_true"),
  pseudo_n_col = c("mean_n_dai_pseudo", "mean_n_dai_pseudo", "mean_n_hani_pseudo", "mean_n_hani_pseudo"),
  rubber_outcome = c("forest_to_cashcrop", "woody_to_cashcrop", "forest_to_cashcrop", "woody_to_cashcrop"),
  farmland_outcome = c("forest_to_farmland", "woody_to_farmland", "forest_to_farmland", "woody_to_farmland")
)

build_panel_rows <- function(map_row) {
  target_coef <- coef_data %>%
    filter(
      source_version == map_row$source_version,
      sample == map_row$sample,
      term == "is_sns"
    )

  target_counts <- count_data %>%
    filter(
      source_version == map_row$source_version,
      sample == map_row$sample,
      stage == "After overlap"
    )

  pseudo_n <- round(target_counts[[map_row$pseudo_n_col]][1])
  true_n <- round(target_counts[[map_row$true_n_col]][1])

  pseudo_row <- tibble(
    panel_id = map_row$panel_id,
    panel_title = map_row$panel_title,
    row_name = "Pseudo",
    outcome = NA_character_,
    estimate = 0,
    ci_low = 0,
    ci_high = 0,
    sig_rate = NA_real_,
    n_value = pseudo_n,
    point_type = "pseudo"
  )

  effect_spec <- tibble(
    row_name = c("To rubber", "To farmland"),
    outcome = c(map_row$rubber_outcome, map_row$farmland_outcome)
  )

  effect_rows <- bind_rows(lapply(seq_len(nrow(effect_spec)), function(i) {
    this_outcome <- effect_spec$outcome[i]
    this_row <- effect_spec$row_name[i]
    draw_df <- target_coef %>%
      filter(outcome == this_outcome)
    tibble(
      panel_id = map_row$panel_id,
      panel_title = map_row$panel_title,
      row_name = this_row,
      outcome = this_outcome,
      estimate = mean(draw_df$estimate, na.rm = TRUE),
      ci_low = quantile(draw_df$estimate, 0.025, na.rm = TRUE),
      ci_high = quantile(draw_df$estimate, 0.975, na.rm = TRUE),
      sig_rate = mean(draw_df$p_value < 0.05, na.rm = TRUE),
      n_value = true_n,
      point_type = "effect"
    )
  }))

  bind_rows(pseudo_row, effect_rows)
}

plot_df <- bind_rows(lapply(split(panel_map, seq_len(nrow(panel_map))), build_panel_rows)) %>%
  mutate(
    row_name = factor(row_name, levels = c("Pseudo", "To rubber", "To farmland")),
    x_pos = c(1, 3, 5)[match(row_name, c("Pseudo", "To rubber", "To farmland"))]
  )

y_min <- -0.30
y_max <- 0.15
global_span <- y_max - y_min
axis_breaks <- c(-0.30, -0.20, -0.10, 0, 0.10, 0.15)

plot_df <- plot_df %>%
  mutate(
    y_n = y_min + global_span * 0.16,
    y_sig = y_max - global_span * 0.04,
    sig_label = ifelse(is.na(sig_rate), "", sprintf("%.2f", sig_rate))
  )

write_csv(plot_df, file.path(out_dir, "repeated_draw_ethnic_effect_summary_data.csv"))

make_single_panel <- function(panel_id_value) {
  panel_df <- plot_df %>%
    filter(panel_id == panel_id_value)

  p <- ggplot(panel_df, aes(x = x_pos, y = estimate)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.8) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "grey75", linewidth = 0.6) +
    geom_segment(
      data = subset(panel_df, point_type == "effect"),
      aes(x = x_pos, xend = x_pos, y = ci_low, yend = ci_high),
      linewidth = 1.0,
      color = "grey55"
    ) +
    geom_point(
      data = subset(panel_df, point_type == "pseudo"),
      shape = 21,
      size = 4.0,
      stroke = 1.0,
      fill = "white",
      color = "grey25"
    ) +
    geom_point(
      data = subset(panel_df, point_type == "effect"),
      shape = 19,
      size = 4.2,
      color = "grey25"
    ) +
    geom_text(
      aes(x = x_pos, y = y_n, label = n_value),
      inherit.aes = FALSE,
      vjust = 1,
      size = 6.2,
      fontface = "bold"
    ) +
    geom_text(
      data = subset(panel_df, point_type == "effect"),
      aes(x = x_pos, y = y_sig, label = sig_label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 6.0,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = c(1, 3, 5),
      labels = rep("", 3),
      limits = c(0.2, 5.8),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = axis_breaks,
      labels = function(x) sprintf("%.0f", x * 100)
    ) +
    labs(
      title = paste0(unique(panel_df$panel_id), "  ", unique(panel_df$panel_title)),
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(ylim = c(y_min, y_max), clip = "on") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 15, color = "black", face = "bold"),
      axis.line.y = element_line(color = "grey35", linewidth = 0.7),
      axis.ticks.y = element_line(color = "grey35", linewidth = 0.7),
      axis.ticks.length.y = unit(5, "pt"),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 24, face = "bold", hjust = 0, margin = margin(b = 8)),
      plot.margin = margin(10, 10, 8, 10)
    )

  file_stub <- paste0("panel_", panel_id_value, "_", gsub(" ", "_", tolower(unique(panel_df$panel_title))))

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".png")),
    plot = p,
    width = 7.2,
    height = 4.8,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".pdf")),
    plot = p,
    width = 7.2,
    height = 4.8,
    bg = "white"
  )
}

invisible(lapply(panel_map$panel_id, make_single_panel))

message("Saved 4 panel-specific plots to: ", out_dir)
