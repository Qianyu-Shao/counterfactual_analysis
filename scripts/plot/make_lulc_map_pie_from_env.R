library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggspatial)
library(scales)
library(purrr)
library(tibble)

# This script assumes these objects already exist in the R environment:
# lulc_2002, lulc_2010, lulc_2014, lulc_2018, lulc_2020, pa

required_objects <- c(
  "lulc_2002", "lulc_2010", "lulc_2014", "lulc_2018", "lulc_2020", "pa"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]

if (length(missing_objects) > 0) {
  stop("Missing objects in environment: ", paste(missing_objects, collapse = ", "))
}

lulc_levels <- c(
  "1" = "Natural forest",
  "2" = "Shrubland",
  "3" = "Rubber plantation",
  "4" = "Tea plantation",
  "5" = "Farmland",
  "6" = "Waterbody",
  "7" = "Others"
)

lulc_name_lookup <- c(
  "natural forests" = "Natural forest",
  "natural forest" = "Natural forest",
  "shrublands" = "Shrubland",
  "shrubland" = "Shrubland",
  "rubber plantations" = "Rubber plantation",
  "rubber plantation" = "Rubber plantation",
  "tea plantations" = "Tea plantation",
  "tea plantation" = "Tea plantation",
  "farmlands" = "Farmland",
  "farmland" = "Farmland",
  "waterbodies" = "Waterbody",
  "waterbody" = "Waterbody",
  "other land uses" = "Others",
  "others" = "Others",
  "other" = "Others"
)

lulc_cols <- c(
  "Natural forest" = "#18A118",
  "Shrubland" = "#F6A313",
  "Rubber plantation" = "#FF1F1F",
  "Tea plantation" = "#FFF100",
  "Farmland" = "#D414F5",
  "Waterbody" = "#52D6F4",
  "Others" = "#000000"
)

output_png <- "./plot/lulc_map_pie_from_env.png"

align_to_template <- function(r, template) {
  if (!same.crs(r, template)) {
    r <- project(r, template, method = "near")
  }

  if (!compareGeom(r, template, stopOnError = FALSE)) {
    r <- resample(r, template, method = "near")
  }

  r
}

raster_to_plot_df <- function(r, year_label, class_lookup) {
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)

  names(df)[3] <- "value"

  if (is.factor(df$value)) {
    df$value <- as.character(df$value)
  }

  if (is.character(df$value)) {
    value_chr <- trimws(df$value)
    suppressWarnings(value_num <- as.integer(value_chr))
    class_from_num <- unname(class_lookup[as.character(value_num)])
    class_from_name <- unname(lulc_name_lookup[tolower(value_chr)])
    df$class <- ifelse(!is.na(class_from_num), class_from_num, class_from_name)
    df$value <- value_num
  } else {
    df$value <- as.integer(round(df$value))
    df$class <- unname(class_lookup[as.character(df$value)])
  }

  df %>%
    mutate(
      class = factor(class, levels = names(lulc_cols)),
      year = year_label
    ) %>%
    filter(!is.na(class)) %>%
    select(x, y, value, class, year)
}

summarise_lulc <- function(df_plot, r, year_label) {
  pixel_area_ha <- prod(res(r)) / 10000

  df_plot %>%
    count(year, class, name = "count") %>%
    mutate(
      area_ha = count * pixel_area_ha,
      prop = area_ha / sum(area_ha)
    )
}

plot_lulc_map <- function(df_plot, pa_sf, year_label, show_north = FALSE, show_scale = FALSE) {
  ggplot() +
    geom_raster(
      data = df_plot,
      aes(x = x, y = y, fill = class)
    ) +
    geom_sf(
      data = pa_sf,
      fill = NA,
      color = "grey35",
      linewidth = 0.35,
      linetype = "22"
    ) +
    scale_fill_manual(
      values = lulc_cols,
      drop = FALSE,
      na.value = "transparent"
    ) +
    coord_sf(expand = FALSE) +
    {
      if (show_north) {
        annotation_north_arrow(
          location = "tl",
          which_north = "true",
          style = north_arrow_fancy_orienteering,
          height = unit(1.0, "cm"),
          width = unit(1.0, "cm")
        )
      }
    } +
    {
      if (show_scale) {
        annotation_scale(location = "bl", width_hint = 0.3)
      }
    } +
    labs(title = year_label, fill = NULL) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    )
}

plot_lulc_pie <- function(df_year) {
  pie_df <- df_year %>%
    mutate(
      class = factor(class, levels = names(lulc_cols))
    ) %>%
    arrange(desc(class)) %>%
    mutate(
      ymax = cumsum(prop),
      ymin = lag(ymax, default = 0),
      ymid = (ymax + ymin) / 2,
      label = percent(prop, accuracy = 0.01),
      hjust = ifelse(ymid < 0.5, 0, 1),
      x_text = ifelse(ymid < 0.5, 2.85, 1.15),
      x_seg_end = ifelse(ymid < 0.5, 2.55, 1.45)
    )

  ggplot(pie_df, aes(x = 2, y = prop, fill = class)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y", clip = "off") +
    scale_fill_manual(values = lulc_cols, drop = FALSE) +
    geom_segment(
      data = subset(pie_df, prop >= 0.02),
      aes(
        x = 2.5,
        xend = x_seg_end,
        y = ymid,
        yend = ymid
      ),
      inherit.aes = FALSE,
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_text(
      data = subset(pie_df, prop >= 0.02),
      aes(
        x = x_text,
        y = ymid,
        label = label,
        hjust = hjust
      ),
      inherit.aes = FALSE,
      size = 4
    ) +
    theme_void() +
    theme(legend.position = "bottom")
}

rasters <- list(
  `2002` = lulc_2002,
  `2010` = lulc_2010,
  `2014` = lulc_2014,
  `2018` = lulc_2018,
  `2020` = lulc_2020
)

template_r <- rasters[[1]]
rasters <- map(rasters, align_to_template, template = template_r)
pa <- st_transform(pa, crs(template_r))

plot_dfs <- imap(rasters, function(r, yr) {
  message("Preparing raster cells for ", yr, "...")
  raster_to_plot_df(r, yr, lulc_levels)
})

lulc_summary <- imap_dfr(rasters, function(r, yr) {
  summarise_lulc(plot_dfs[[yr]], r, yr)
})

years <- names(rasters)

map_plots <- map(
  years,
  function(yr) {
    plot_lulc_map(
      df_plot = plot_dfs[[yr]],
      pa_sf = pa,
      year_label = yr,
      show_north = (yr == years[1]),
      show_scale = (yr == tail(years, 1))
    )
  }
)

pie_plots <- map(years, ~ plot_lulc_pie(filter(lulc_summary, year == .x)))

row_plots <- map2(
  map_plots,
  pie_plots,
  ~ .x + .y + plot_layout(widths = c(1.3, 0.9))
)

final_plot <- wrap_plots(row_plots, ncol = 1)
final_plot <- final_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(final_plot)

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 12,
  height = 4 * length(years),
  dpi = 300,
  bg = "white"
)

message("Saved figure to: ", output_png)
