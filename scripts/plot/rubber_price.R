

#### general distribution of sacred sites of different ethnic groups

sacred_points <- sacred_points %>% mutate(ethnic_group = case_when(
  民族 %in% c("水傣", "汉傣", "花腰傣") & 有无寺庙 == "有" ~ "Dai with Temple",
  民族 %in% c("水傣", "汉傣", "花腰傣") & 有无寺庙 == "无" ~ "Dai without Temple",
  民族 %in% c("哈尼") ~ "Hani",
  民族 %in% c("布朗") ~ "Bulang"
)) %>% 
  filter(!is.na(ethnic_group))


library(sf)
library(terra)
library(ggplot2)
library(tidyterra) 
library(ggnewscale)

# 1. 加载数据
# 假设 srtm_raster 是你的 SRTM SpatRaster，sns_pts 是 930 个点的 sf 对象
# srtm <- rast("path_to_srtm.tif")
# sns_pts <- st_read("path_to_points.shp")

# 2. 准备底图：计算坡向和坡度以生成阴影图 (可选，但非常美观)
slope <- terrain(srtm, "slope", unit = "radians")
aspect <- terrain(srtm, "aspect", unit = "radians")
hill <- shade(slope, aspect, angle = 45, direction = 315)

# 3. 绘图
ggplot() +
  # 绘制山体阴影作为底层
  geom_spatraster(data = hill, show.legend = FALSE) +
  scale_fill_gradient(low = "grey30", high = "white") +
  
  # 叠加 SRTM 高程（设置半透明，让阴影透出来）
  new_scale_fill() + # 开启新的 fill 标尺，防止冲突
  geom_spatraster(data = srtm, alpha = 0.4) +
  scale_fill_terrain_c(name = "Elevation (m)") +
  
  # 叠加圣境点，按民族分组
  geom_sf(data = sacred_points, aes(color = ethnic_group), size = 2, alpha = 0.8) +
  
  # 设置颜色方案
  scale_color_brewer(palette = "Set1", name = "Ethnic Group") +
  
  # 细节修饰
  theme_minimal() +
  labs(
    title = "Spatial Distribution of Sacred Sites by Ethnic Group",
    subtitle = "Overlay on SRTM Terrain",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(legend.position = "right")




full_path <- "./plot/distribution_of_sites.png"

ggsave(
  filename = full_path,
  plot = last_plot(),  # 保存最后一次显示的图
  width = 10,          # 宽度（默认单位为英寸）
  height = 8,          # 高度
  dpi = 300,           # 分辨率，300是学术论文的标准
  bg = "white"         # 背景色，防止透明背景导致黑图
)




###################rubber price
library(rvest)
library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

file_path <- "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis/data/rubber price/rubber-300.xls"

tbl_list <- html_table(read_html(file_path))
rubber_raw <- tbl_list[[1]]

rubber_monthly <- rubber_raw %>%
  transmute(
    month = my(Month),
    price = as.numeric(Price)
  ) %>%
  filter(!is.na(month), !is.na(price)) %>%
  filter(year(month) >= 2002, year(month) <= 2020) %>%
  arrange(month)

p <- ggplot(rubber_monthly, aes(x = month, y = price)) +
  geom_line(linewidth = 0.7, color = "black") +
  scale_x_date(
    limits = as.Date(c("2002-01-01", "2020-12-31")),
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Price (USD/kg)",
    title = "Monthly Natural Rubber Price, 2002-2020"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

p



ggsave(
  "./plot/rubber_price_2002_2020.png",
  plot = p,
  width = 8,
  height = 4.5,
  dpi = 300,
  bg = "white"
)


