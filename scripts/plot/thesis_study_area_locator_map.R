suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(tidyterra)
  library(ggnewscale)
  library(patchwork)
  library(ggspatial)
  library(grid)
})

project_root <- "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis"
setwd(project_root)

invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

china_lvl0 <- st_read("./data/gadm41_CHN_shp/gadm41_CHN_0.shp", quiet = TRUE)
china_lvl1 <- st_read("./data/gadm41_CHN_shp/gadm41_CHN_1.shp", quiet = TRUE)
china_lvl2 <- st_read("./data/gadm41_CHN_shp/gadm41_CHN_2.shp", quiet = TRUE)
china_border <- st_read("./data/china-geospatial-data-UTF8/CN-border-La.gmt", quiet = TRUE)
china_border_lines <- st_read("./data/china-geospatial-data-UTF8/CN-border-L1.gmt", quiet = TRUE)
ten_dash_line <- st_read("./data/china-geospatial-data-UTF8/ten-dash-line.gmt", quiet = TRUE) %>%
  st_transform(4326)
xsbn_prefecture_boundary <- st_read(
  "./data/vector/Xishuangbanna_Boundaries/Prefecture_boundary.shp",
  quiet = TRUE
) %>%
  st_transform(4326)

yunnan <- china_lvl1 %>%
  filter(NAME_1 == "Yunnan")

yunnan_locator <- china_border %>%
  filter(CODE == 530000)

china_provinces <- china_border %>%
  filter(!is.na(CODE))

world_map <- map_data("world") %>%
  filter(
    long >= 70,
    long <= 138
  )

xsbn <- china_lvl2 %>%
  filter(NAME_1 == "Yunnan", grepl("Xishuangbanna", NAME_2, ignore.case = TRUE))

sacred_points <- st_read(
  "./data/vector/sacred/自然圣境数据_西双版纳-yuanshi.shp",
  quiet = TRUE
) %>%
  st_transform(4326) %>%
  mutate(
    ethnic_group = case_when(
      民族 %in% c("水傣", "汉傣", "花腰傣") ~ "Dai",
      民族 == "哈尼" ~ "Hani",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ethnic_group)) %>%
  mutate(
    ethnic_group = factor(ethnic_group, levels = c("Dai", "Hani"))
  )

srtm <- rast("./data/raster/banna_SRTM.tif")
slope <- terrain(srtm, "slope", unit = "radians")
aspect <- terrain(srtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, angle = 45, direction = 315)

study_bbox <- st_as_sfc(st_bbox(srtm))
study_bbox_sf <- st_as_sf(study_bbox)
st_crs(study_bbox_sf) <- st_crs(srtm)
study_bbox_ll <- st_transform(study_bbox_sf, 4326)

xsbn_bbox <- st_bbox(xsbn_prefecture_boundary)
xsbn_xlim <- c(xsbn_bbox["xmin"] - 0.08, xsbn_bbox["xmax"] + 0.08)
xsbn_ylim <- c(xsbn_bbox["ymin"] - 0.08, xsbn_bbox["ymax"] + 0.08)

make_hatch_lines <- function(poly, spacing = 0.13) {
  bbox <- st_bbox(poly)
  width <- as.numeric(bbox["xmax"] - bbox["xmin"])
  starts <- seq(
    from = as.numeric(bbox["ymin"] - width),
    to = as.numeric(bbox["ymax"]),
    by = spacing
  )
  lines <- st_sfc(
    lapply(starts, function(y0) {
      st_linestring(matrix(
        c(
          bbox["xmin"], y0,
          bbox["xmax"], y0 + width
        ),
        ncol = 2,
        byrow = TRUE
      ))
    }),
    crs = st_crs(poly)
  )
  suppressWarnings(st_intersection(st_sf(geometry = lines), st_union(poly)))
}

yunnan_hatch <- make_hatch_lines(yunnan_locator, spacing = 0.70)
xsbn_vect <- vect(xsbn_prefecture_boundary)
srtm_masked <- mask(crop(srtm, xsbn_vect), xsbn_vect)
hill_masked <- mask(crop(hill, xsbn_vect), xsbn_vect)

elevation_palette <- c(
  "#2F7D32", "#78B943", "#D6DD57", "#F2CF72", "#B77A45", "#E6D8D0"
)

format_lon <- function(x) paste0(round(x), "\u00b0E")
format_lat <- function(x) {
  paste0(abs(round(x)), "\u00b0", ifelse(x < 0, "S", "N"))
}

locator_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

p_china <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "#F3F0E6",
    color = "#D8D2C4",
    linewidth = 0.10
  ) +
  geom_sf(data = china_border, fill = "#FAF7EF", color = NA) +
  geom_sf(data = china_provinces, fill = NA, color = "#D9D3C7", linewidth = 0.13) +
  geom_sf(data = china_border_lines, color = "black", linewidth = 0.45) +
  geom_sf(
    data = yunnan_locator,
    fill = NA,
    color = "#7A5A48",
    linewidth = 0.30
  ) +
  geom_sf(
    data = xsbn_prefecture_boundary,
    fill = "#9CC9C7",
    color = "#4F8F8C",
    linewidth = 0.24
  ) +
  geom_sf(data = yunnan_hatch, color = "#7A5A48", linewidth = 0.18) +
  geom_sf(data = ten_dash_line, color = "black", linewidth = 0.50) +
  annotate(
    "text",
    x = 104,
    y = 37,
    label = "China",
    color = "grey25",
    fontface = "bold",
    size = 4.2
  ) +
  annotation_scale(
    location = "bl",
    style = "ticks",
    width_hint = 0.20,
    text_cex = 0.65,
    line_width = 0.45,
    pad_x = unit(0.42, "cm"),
    pad_y = unit(0.35, "cm")
  ) +
  annotate(
    "rect",
    xmin = 78.5,
    xmax = 81.3,
    ymin = 12.4,
    ymax = 15.0,
    fill = NA,
    color = "#7A5A48",
    linewidth = 0.25
  ) +
  annotate(
    "segment",
    x = c(78.6, 79.4, 80.2, 81.0),
    xend = c(79.5, 80.3, 81.1, 81.3),
    y = c(12.5, 12.5, 12.5, 13.1),
    yend = c(14.9, 14.9, 14.9, 14.9),
    color = "#7A5A48",
    linewidth = 0.25
  ) +
  annotate(
    "text",
    x = 82.6,
    y = 13.7,
    label = "Yunnan Province",
    hjust = 0,
    size = 3.4
  ) +
  annotate(
    "rect",
    xmin = 78.5,
    xmax = 81.3,
    ymin = 8.6,
    ymax = 11.2,
    fill = "#9CC9C7",
    color = "#4F8F8C",
    linewidth = 0.25
  ) +
  annotate(
    "text",
    x = 82.6,
    y = 9.9,
    label = "Xishuangbanna",
    hjust = 0,
    size = 3.4
  ) +
  coord_sf(
    xlim = c(73, 135),
    ylim = c(-12, 64),
    expand = FALSE
  ) +
  scale_x_continuous(
    position = "top",
    breaks = c(80, 100, 120),
    labels = format_lon
  ) +
  scale_y_continuous(
    position = "left",
    breaks = c(-10, 10, 30, 50),
    labels = format_lat
  ) +
  labs(title = NULL) +
  locator_theme(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

ggsave(
  filename = "./plot/china_map_yunnan_highlight.png",
  plot = p_china,
  width = 5.6,
  height = 6.8,
  dpi = 500,
  bg = "white"
)

p_study <- ggplot() +
  geom_spatraster(data = hill_masked, show.legend = FALSE) +
  scale_fill_gradient(low = "grey30", high = "white", na.value = NA) +
  new_scale_fill() +
  geom_spatraster(data = srtm_masked, alpha = 0.70) +
  scale_fill_gradientn(
    name = "Elevation (m)",
    colors = elevation_palette,
    limits = c(0, 2420),
    breaks = c(0, 1200, 2400),
    labels = c("0", "1,200", "2,400"),
    na.value = NA,
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(3.8, "cm"),
      barheight = unit(0.30, "cm"),
      ticks.colour = "grey25"
    )
  ) +
  geom_sf(
    data = sacred_points,
    aes(color = ethnic_group),
    size = 1.75,
    alpha = 0.90
  ) +
  geom_sf(
    data = xsbn_prefecture_boundary,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  annotation_scale(
    location = "bl",
    style = "ticks",
    width_hint = 0.18,
    text_cex = 0.7,
    line_width = 0.5,
    pad_x = unit(0.55, "cm"),
    pad_y = unit(0.35, "cm")
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      line_col = "black",
      text_col = "black"
    ),
    height = unit(0.85, "cm"),
    width = unit(0.85, "cm"),
    pad_x = unit(0.35, "cm"),
    pad_y = unit(0.35, "cm")
  ) +
  scale_color_manual(
    values = c("Dai" = "#2166AC", "Hani" = "#7B3294"),
    name = "Ethnic Group",
    guide = guide_legend(
      title.position = "top",
      nrow = 1,
      override.aes = list(size = 2.6, alpha = 1)
    )
  ) +
  coord_sf(
    xlim = xsbn_xlim,
    ylim = xsbn_ylim,
    expand = FALSE,
    label_axes = list(bottom = "E", right = "N")
  ) +
  scale_x_continuous(
    position = "bottom",
    breaks = c(100.0, 100.5, 101.0, 101.5),
    labels = function(x) paste0(x, "\u00b0E")
  ) +
  scale_y_continuous(
    position = "right",
    breaks = c(21.5, 22.0, 22.5),
    labels = function(x) paste0(x, "\u00b0N")
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.text = element_text(color = "black", face = "bold", size = 9),
    legend.position = c(0.33, 0.12),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.background = element_rect(fill = scales::alpha("white", 0.72), color = NA),
    legend.key = element_rect(fill = NA, color = NA),
    legend.spacing.x = unit(0.22, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

final_plot <- (free(p_china) | free(p_study)) + plot_layout(widths = c(0.92, 1.48))

final_plot <- final_plot +
  plot_annotation(
    title = "Location of the Study Area and Spatial Distribution of Sacred Sites",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )

output_path <- "./plot/thesis_study_area_locator_map.png"

ggsave(
  filename = output_path,
  plot = final_plot,
  width = 12.8,
  height = 7.6,
  dpi = 400,
  bg = "white"
)

message("Saved:")
message(output_path)
