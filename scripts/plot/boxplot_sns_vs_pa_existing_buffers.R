`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MatchIt)
  library(rstatix)
  library(ggpubr)
})

if (!exists("final_buffer_500_12345_forest") || !exists("final_buffer_500_12345_woody")) {
  stop(
    "final_buffer_500_12345_forest and final_buffer_500_12345_woody must already exist in the current R session."
  )
}

source("./scripts/function.R")

if (!exists("lulc_2002")) {
  lulc_2002 <- terra::rast("./data/raster/LULC/XSBN_LULC_02_14/2002-lucc-ff.img")
}
if (!exists("lulc_2014")) {
  lulc_2014 <- terra::rast("./data/raster/LULC/XSBN_LULC_02_14/xsbn-2014-final_recl.img")
}

target_crs <- "EPSG:32647"
lulc_2002 <- terra::project(lulc_2002, target_crs, method = "near")
lulc_2014 <- terra::project(lulc_2014, target_crs, method = "near")

dir.create("./plot", showWarnings = FALSE, recursive = TRUE)

output_plot <- "./plot/sns_vs_pa_existing_buffers_within_ethnicity_boxplot.png"
output_counts <- "./plot/sns_vs_pa_existing_buffers_within_ethnicity_counts.csv"

ensure_sf_x <- function(obj) {
  if (inherits(obj, "sf")) {
    if (attr(obj, "sf_column") != "x" && "x" %in% names(obj)) {
      sf::st_geometry(obj) <- "x"
    }
    return(obj)
  }

  if (!"x" %in% names(obj)) {
    stop("Object must contain a geometry column named `x`.")
  }

  sf::st_as_sf(obj, sf_column_name = "x", crs = target_crs)
}

forest_sf <- ensure_sf_x(final_buffer_500_12345_forest)
woody_sf <- ensure_sf_x(final_buffer_500_12345_woody)

prepare_source <- function(sf_data) {
  sf_data %>%
    mutate(
      ethnicity = case_when(
        is_sns == 1 & dai == 1 ~ "Dai",
        is_sns == 1 & hani == 1 ~ "Hani",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(is_sns == 0 | ethnicity %in% c("Dai", "Hani"))
}

forest_source <- prepare_source(forest_sf)
woody_source <- prepare_source(woody_sf)

run_match_safe <- function(sf_data, eth_name, primary_caliper = 0.05) {
  calipers <- unique(c(primary_caliper, 0.05, 0.08, 0.10, 0.15))

  match_df <- sf_data %>%
    st_drop_geometry() %>%
    mutate(
      ethnicity = replace_na(ethnicity, "PA"),
      aspect_sin = sin(mode_aspect * pi / 180),
      aspect_cos = cos(mode_aspect * pi / 180)
    ) %>%
    filter(is_sns == 0 | ethnicity == eth_name)

  for (cal in calipers) {
    match_obj <- tryCatch(
      matchit(
        is_sns ~ mean_elev + mean_slope + aspect_sin + aspect_cos,
        data = match_df,
        method = "nearest",
        distance = "glm",
        caliper = cal,
        replace = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(match_obj)) {
      next
    }

    matched_df <- tryCatch(match.data(match_obj), error = function(e) NULL)
    if (is.null(matched_df) || nrow(matched_df) == 0) {
      next
    }

    matched_sf <- sf_data %>%
      inner_join(
        matched_df %>%
          select(ID, is_sns, distance, weights, subclass),
        by = c("ID", "is_sns")
      ) %>%
      mutate(ethnicity = eth_name)

    return(list(
      match_obj = match_obj,
      matched_data = matched_sf,
      used_caliper = cal
    ))
  }

  stop(paste("No units were matched for", eth_name))
}

forest_match_dai <- run_match_safe(forest_source, "Dai", primary_caliper = 0.05)
forest_match_hani <- run_match_safe(forest_source, "Hani", primary_caliper = 0.05)

woody_match_dai <- run_match_safe(woody_source, "Dai", primary_caliper = 0.05)
woody_match_hani <- run_match_safe(woody_source, "Hani", primary_caliper = 0.05)

matched_forest_pa_dai_hani <- bind_rows(
  forest_match_dai$matched_data,
  forest_match_hani$matched_data
)

 matched_forest_pa_dai <- forest_match_dai$matched_data


 st_write(matched_forest_pa_dai, "matched_forest_pa_dai.gpkg", delete_dsn = TRUE)


matched_woody_pa_dai_hani <- bind_rows(
  woody_match_dai$matched_data,
  woody_match_hani$matched_data
)

# Add 2002-2014 intensity after matching
matched_forest_pa_dai_hani$loss_1_34_14 <- calculate_lulc_loss_1_34(
  lulc_2002, lulc_2014, matched_forest_pa_dai_hani
)
matched_forest_pa_dai_hani$loss_1_34_intensity_14 <-
  matched_forest_pa_dai_hani$loss_1_34_14 / matched_forest_pa_dai_hani$forest_area_02

matched_forest_pa_dai_hani$loss_1_345_14 <- calculate_lulc_loss_1_5(
  lulc_2002, lulc_2014, matched_forest_pa_dai_hani
)
matched_forest_pa_dai_hani$loss_1_345_intensity_14 <-
  matched_forest_pa_dai_hani$loss_1_345_14 / matched_forest_pa_dai_hani$forest_area_02

matched_woody_pa_dai_hani$loss_12_34_14 <- calculate_lulc_loss_12_34(
  lulc_2002, lulc_2014, matched_woody_pa_dai_hani
)
matched_woody_pa_dai_hani$loss_12_34_intensity_14 <-
  matched_woody_pa_dai_hani$loss_12_34_14 / matched_woody_pa_dai_hani$woody_area_02

matched_woody_pa_dai_hani$loss_12_345_14 <- calculate_lulc_loss_12_5(
  lulc_2002, lulc_2014, matched_woody_pa_dai_hani
)
matched_woody_pa_dai_hani$loss_12_345_intensity_14 <-
  matched_woody_pa_dai_hani$loss_12_345_14 / matched_woody_pa_dai_hani$woody_area_02

assign("matched_forest_pa_dai_hani", matched_forest_pa_dai_hani, envir = .GlobalEnv)
assign("matched_woody_pa_dai_hani", matched_woody_pa_dai_hani, envir = .GlobalEnv)

extract_counts <- function(source_df, matched_df, sample_name, eth_name) {
  before <- source_df %>%
    st_drop_geometry() %>%
    filter(is_sns == 0 | ethnicity == eth_name) %>%
    mutate(group = ifelse(is_sns == 1, "True SNS", "PA pseudo")) %>%
    count(group, name = "n") %>%
    mutate(stage = "Before")

  after <- matched_df %>%
    st_drop_geometry() %>%
    mutate(group = ifelse(is_sns == 1, "True SNS", "PA pseudo")) %>%
    count(group, name = "n") %>%
    mutate(stage = "After")

  bind_rows(before, after) %>%
    mutate(sample = sample_name, ethnicity = eth_name) %>%
    select(sample, ethnicity, stage, group, n)
}

count_table <- bind_rows(
  extract_counts(forest_source, forest_match_dai$matched_data, "Forest", "Dai"),
  extract_counts(forest_source, forest_match_hani$matched_data, "Forest", "Hani"),
  extract_counts(woody_source, woody_match_dai$matched_data, "Shrubland", "Dai"),
  extract_counts(woody_source, woody_match_hani$matched_data, "Shrubland", "Hani")
)

write.csv(count_table, output_counts, row.names = FALSE)

make_ethnic_long <- function(forest_df, woody_df, eth_name) {
  forest_plot <- forest_df %>%
    filter(ethnicity == eth_name) %>%
    st_drop_geometry() %>%
    select(is_sns, loss_1_34_intensity_14, loss_1_345_intensity_14) %>%
    pivot_longer(
      cols = c(loss_1_34_intensity_14, loss_1_345_intensity_14),
      names_to = "outcome",
      values_to = "encroachment"
    )

  woody_plot <- woody_df %>%
    filter(ethnicity == eth_name) %>%
    st_drop_geometry() %>%
    select(is_sns, loss_12_34_intensity_14, loss_12_345_intensity_14) %>%
    pivot_longer(
      cols = c(loss_12_34_intensity_14, loss_12_345_intensity_14),
      names_to = "outcome",
      values_to = "encroachment"
    )

  bind_rows(forest_plot, woody_plot) %>%
    mutate(ethnicity = eth_name)
}

plot_data <- bind_rows(
  make_ethnic_long(matched_forest_pa_dai_hani, matched_woody_pa_dai_hani, "Dai"),
  make_ethnic_long(matched_forest_pa_dai_hani, matched_woody_pa_dai_hani, "Hani")
) %>%
  mutate(
    ethnicity = factor(ethnicity, levels = c("Dai", "Hani")),
    is_sns = factor(is_sns, levels = c(0, 1), labels = c("PA pseudo", "True SNS")),
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
  filter(!is.na(encroachment))

stat_test <- plot_data %>%
  group_by(ethnicity, outcome) %>%
  wilcox_test(encroachment ~ is_sns) %>%
  ungroup()

direction_df <- plot_data %>%
  group_by(ethnicity, outcome, is_sns) %>%
  summarise(mean_enc = mean(encroachment, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = is_sns, values_from = mean_enc) %>%
  mutate(
    diff = `True SNS` - `PA pseudo`,
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
  scale_fill_manual(values = c("PA pseudo" = "#D95F02", "True SNS" = "#1B9E77")) +
  scale_color_manual(values = c("PA pseudo" = "#D95F02", "True SNS" = "#1B9E77")) +
  labs(
    title = "Encroachment Intensity After Ethnicity-Specific Matching of SNS and PA Sites",
    subtitle = "\u2191 True SNS > PA pseudo; \u2193 True SNS < PA pseudo",
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

ggsave(output_plot, p, width = 16, height = 10, dpi = 300)

message("Saved plot: ", output_plot)
message("Saved counts: ", output_counts)



summary(forest_match_dai$match_obj)
summary(forest_match_hani$match_obj)

summary(woody_match_dai$match_obj)
summary(woody_match_hani$match_obj)
