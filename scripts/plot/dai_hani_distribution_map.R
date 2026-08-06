suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(tidyterra)
  library(ggnewscale)
})

project_root <- "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis"
setwd(project_root)

invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

sacred_points <- st_read(
  "./data/vector/sacred/自然圣境数据_西双版纳-yuanshi.shp",
  quiet = TRUE
) %>%
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

p <- ggplot() +
  geom_spatraster(data = hill, show.legend = FALSE) +
  scale_fill_gradient(low = "grey30", high = "white") +
  new_scale_fill() +
  geom_spatraster(data = srtm, alpha = 0.42) +
  scale_fill_terrain_c(name = "Elevation (m)") +
  geom_sf(
    data = sacred_points,
    aes(color = ethnic_group),
    size = 2.1,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c("Dai" = "#2C7FB8", "Hani" = "#984EA3"),
    name = "Ethnic Group"
  ) +
  labs(
    title = "Spatial Distribution of Dai and Hani Sacred Sites",
    subtitle = "Overlay on SRTM Terrain",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

output_path <- "./plot/dai_hani_distribution_of_sites.png"

ggsave(
  filename = output_path,
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("Saved:")
message(output_path)
