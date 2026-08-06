
library(sf)
library(terra)
library(dplyr)
library(tidyr)

library(ggplot2)
library(ggspatial)
library(scales)

library(readr)
library(readxl)
library(units)

library(MatchIt)
library(spdep)
library(spatialreg)

library(car)



# 1. 读取数据 


village <- st_read("./data/vector/village/village.shp")

sacred_points <- st_read("./data/vector/sacred/自然圣境数据_西双版纳-yuanshi.shp")

lulc_2002 <- rast("./data/raster/LULC/XSBN_LULC_02_14/2002-lucc-ff.img")


lulc_2014 <- rast("./data/raster/LULC/XSBN_LULC_02_14/xsbn-2014-final_recl.img")

lulc_2020 <- rast("./data/raster/LULC/lulc_2020/lc_2020_xsbn_reclas.tif")

srtm <- rast("./data/raster/banna_SRTM.tif")

social_village <- read_excel("./data/eco_data.xlsx", sheet = "village1314_income_adjusted")


pa <- st_read("./data/vector/pa/XSBN_PAs.shp")

road <- st_read("./data/vector/road/road.shp")

village_pts <- st_read("./data/vector/YN-villages/yn_villages.shp")

# 2. 确保投影


target_crs <- "EPSG:32647"

sacred_points <- st_transform(sacred_points, target_crs)
pa <- st_transform(pa, target_crs)
road <- st_transform(road, target_crs)
village <- st_transform(village, target_crs)
village_pts <- st_transform(village_pts, target_crs)

srtm <- project(srtm, target_crs)
lulc_2002 <- project(lulc_2002, target_crs, method = "near")
lulc_2014 <- project(lulc_2014, target_crs, method = "near")

# 3. 合并数据、删除na


sacred_points$ID <- 1:nrow(sacred_points)
sacred_points <- sacred_points %>% select(ID, everything())




# 计算村庄面积

village$area_sqm <- st_area(village)

village$area_km2 <- drop_units(village$area_sqm) / 1000000



village_clean_1 <- village %>%
  left_join(social_village %>% 
              st_drop_geometry() %>% 
              dplyr::select(OBJECTID, net_income_adjusted, popu, agri_land, dai, hani), 
            by = "OBJECTID") %>%
  filter(!is.na(net_income_adjusted), !is.na(popu), !is.na(agri_land))

village_clean <- village %>%
  left_join(social_village %>% 
              st_drop_geometry() %>% 
              dplyr::select(OBJECTID, net_income_adjusted, popu, agri_land, female, dai, hani, bulang), 
            by = "OBJECTID") %>%
  filter(!is.na(net_income_adjusted), !is.na(popu), !is.na(agri_land), !is.na(dai), !is.na(hani), !is.na(bulang))



village_clean <- village_clean %>%
  mutate(popu_density = popu / area_km2,
         agri_land_pp = agri_land / popu) 



village_clean_1 <- village_clean_1 %>%
  mutate(popu_density = popu / area_km2,
         agri_land_pp = agri_land / popu) 










# 4. 生成buffer

# 筛选：仅保留与村庄边界有交集（即落在里面）的点
pts_in_village <- sacred_points
buffers_500 <- st_buffer(pts_in_village, dist = 500)

# 移除位于pa中的点
intersect_logical_500 <- st_intersects(buffers_500, st_union(pa), sparse = FALSE)[, 1]
final_buffers_true_500 <- buffers_500[!intersect_logical_500, ]

final_buffers_true_500 <- final_buffers_true_500 %>%
  mutate(is_sns = 1)

# st_write(final_buffers, "final_buffer_500m.gpkg", delete_dsn = TRUE)
# st_write(village_clean, "village_clean.gpkg", delete_dsn = TRUE)





# 在村庄创建pseudo buffer (且lulc_2002 = 1/2/3/4/5)

buffers_1000m <- st_buffer(sacred_points, dist = 1000)




final_10000_pts_12345_500 <- select_stratified_points(
  study_area = pa,
  lulc_raster = lulc_2002,
  target_lulc = c(1, 2, 3, 4, 5),  # 在这里更改点下的 LULC 类型
  exclude_polygons_list = list(buffers_1000m), # 放入所有要剔除的层
  initial_size = 30000,
  final_size = 10000
)


# st_write(final_10000_pts_12345_500, "final_10000_pts_12345_500.gpkg", delete_dsn = TRUE)








# 5. 提取海拔、朝向、坡度、道路距离

sacred_true_points_500 <- sacred_points %>%
  right_join(final_buffers_true_500 %>% 
               st_drop_geometry() %>% 
               dplyr::select(ID), 
             by = "ID") %>%
  mutate(is_sns = 1) %>%
  select("is_sns", "民族", "丧葬类型", "有无寺庙") %>%
  rename(x = geometry,
         ethnicity = 民族,
         burial = 丧葬类型,
         temple = 有无寺庙) %>%
  mutate(dai = ifelse(ethnicity %in% c("水傣", "汉傣", "花腰傣"), 1, 0),
         hani = ifelse(ethnicity == "哈尼", 1, 0),
         bulang = ifelse(ethnicity == "布朗", 1, 0))


sacred_false_points_12345_500 <- final_10000_pts_12345_500 %>%
  mutate(is_sns = 0)



total_data_points_12345_500 <- bind_rows(sacred_true_points_500, sacred_false_points_12345_500)



total_data_points_12345_500$ID <- 1:nrow(total_data_points_12345_500)



# 提取公路距离
road_union <- st_union(road)

dist_matrix_12345_500 <- st_distance(total_data_points_12345_500, road_union)

total_data_points_12345_500$dist_to_road <- as.numeric(dist_matrix_12345_500)


# 提取离最近村庄数据


total_data_points_12345_500$dist_to_village <- get_dist_vector(total_data_points_12345_500, village_pts, village)




# 提取海拔、朝向、坡度



total_buffer_500_12345 <- st_buffer(total_data_points_12345_500, dist = 500)

final_buffer_500_12345 <- extract_buffer_terrain(srtm, total_buffer_500_12345)




# 6. 计算initial_forest/woody land 和 encroachment

final_buffer_500_12345$forest_area_02 <- calculate_lulc_area(lulc_2002, final_buffer_500_12345, c(1))
final_buffer_500_12345$woody_area_02 <- calculate_lulc_area(lulc_2002, final_buffer_500_12345, c(2))


final_buffer_500_12345_forest <- final_buffer_500_12345 %>%
  filter(forest_area_02 > 0)

final_buffer_500_12345_woody <- final_buffer_500_12345 %>%
  filter(woody_area_02 > 0)






