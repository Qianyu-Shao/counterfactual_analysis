# ============================================================
# Geography Profile of All Dai and Hani True Sacred Sites
# Purpose:
# 1. Read original sacred points directly
# 2. Extract terrain information from SRTM using 500 m buffers
# 3. Visualize Dai and Hani true sacred-site geography distributions
#    using the same style as Part 1 of
#    dai_hani_geography_vs_institution_analysis.R
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
})

invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

source("./scripts/function.R")

target_crs <- "EPSG:32647"

sacred_points <- st_read(
  "./data/vector/sacred/自然圣境数据_西双版纳-yuanshi.shp"
) %>%
  st_transform(target_crs)

srtm <- rast("./data/raster/banna_SRTM.tif") %>%
  project(target_crs)

make_aspect_sector <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x >= 337.5 | x < 22.5 ~ "North",
    x >= 22.5 & x < 67.5 ~ "Northeast",
    x >= 67.5 & x < 112.5 ~ "East",
    x >= 112.5 & x < 157.5 ~ "Southeast",
    x >= 157.5 & x < 202.5 ~ "South",
    x >= 202.5 & x < 247.5 ~ "Southwest",
    x >= 247.5 & x < 292.5 ~ "West",
    x >= 292.5 & x < 337.5 ~ "Northwest"
  )
}

true_dh_points <- sacred_points %>%
  mutate(
    ethnicity = 民族,
    dai = ifelse(ethnicity %in% c("水傣", "汉傣", "花腰傣"), 1, 0),
    hani = ifelse(ethnicity == "哈尼", 1, 0)
  ) %>%
  filter(dai == 1 | hani == 1) %>%
  mutate(
    ethnicity = ifelse(dai == 1, "Dai", "Hani")
  )

true_dh_buffer_500 <- st_buffer(true_dh_points, dist = 500)

true_dh_terrain_500 <- extract_buffer_terrain(
  srtm = srtm,
  buffer_vec = true_dh_buffer_500
) %>%
  mutate(
    ethnicity = factor(ethnicity, levels = c("Dai", "Hani")),
    aspect_sector = make_aspect_sector(mode_aspect)
  )

ethnicity_colors <- c("Dai" = "#D4D4CC", "Hani" = "#86AE8E")
ethnicity_line_colors <- c("Dai" = "#85857E", "Hani" = "#3F7A4D")

plot_density_with_points <- function(df, file_name) {
  density_df <- df %>%
    st_drop_geometry() %>%
    select(ethnicity, mean_elev, mean_slope) %>%
    pivot_longer(
      cols = c(mean_elev, mean_slope),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      variable = dplyr::recode(
        variable,
        mean_elev = "Elevation (m)",
        mean_slope = "Slope (deg)"
      )
    )

  p <- ggplot(
    density_df,
    aes(x = value, color = ethnicity, fill = ethnicity)
  ) +
    geom_density(alpha = 0.12, adjust = 1.1, linewidth = 0.9) +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    scale_color_manual(values = ethnicity_line_colors) +
    scale_fill_manual(values = ethnicity_colors) +
    labs(
      title = "Geographic Distributions of Dai and Hani True Sacred Sites",
      x = NULL,
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  ggsave(file_name, p, width = 10, height = 5.5, dpi = 300)
  p
}

plot_aspect <- function(df, file_name) {
  aspect_df <- df %>%
    st_drop_geometry() %>%
    count(ethnicity, aspect_sector) %>%
    filter(!is.na(aspect_sector)) %>%
    mutate(
      aspect_sector = factor(
        aspect_sector,
        levels = c(
          "North", "Northeast", "East", "Southeast",
          "South", "Southwest", "West", "Northwest"
        )
      )
    )

  p <- ggplot(aspect_df, aes(x = aspect_sector, y = n, fill = ethnicity)) +
    geom_col(alpha = 0.95) +
    coord_polar() +
    facet_wrap(~ ethnicity, nrow = 1) +
    scale_fill_manual(values = ethnicity_colors) +
    labs(
      title = "Aspect Distribution of Dai and Hani True Sacred Sites",
      x = NULL,
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(file_name, p, width = 11, height = 6, dpi = 300)
  p
}

p_density <- plot_density_with_points(
  true_dh_terrain_500,
  "./plot/all_dai_hani_true_sites_density.png"
)

p_aspect <- plot_aspect(
  true_dh_terrain_500,
  "./plot/all_dai_hani_true_sites_aspect.png"
)

true_geo_compare <- true_dh_terrain_500 %>%
  st_drop_geometry() %>%
  select(ethnicity, mean_elev, mean_slope) %>%
  pivot_longer(
    cols = c(mean_elev, mean_slope),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  wilcox_test(value ~ ethnicity) %>%
  ungroup()

write.csv(
  true_geo_compare,
  "./plot/all_dai_hani_true_sites_geography_compare.csv",
  row.names = FALSE
)

message("Saved:")
message("./plot/all_dai_hani_true_sites_density.png")
message("./plot/all_dai_hani_true_sites_aspect.png")
message("./plot/all_dai_hani_true_sites_geography_compare.csv")
