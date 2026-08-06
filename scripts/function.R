
library(terra)
library(sf)
library(exactextractr)
library(dplyr)





extract_buffer_terrain <- function(srtm, buffer_vec) {
  
  # 1. 统一坐标系 (使用 Buffer 的投影，通常是 UTM)
  target_crs <- crs(buffer_vec)
  if (crs(srtm) != target_crs) {
    message("正在重投影 SRTM 以匹配缓冲区坐标系...")
    srtm <- project(srtm, target_crs)
  }
  
  # 2. 生成地形衍生图层
  message("正在生成坡度与坡向图层...")
  slope <- terrain(srtm, v = "slope", unit = "degrees")
  aspect <- terrain(srtm, v = "aspect", unit = "degrees")
  
  # 合并为栅格堆栈
  terrain_stack <- c(srtm, slope, aspect)
  names(terrain_stack) <- c("elev", "slope", "aspect")
  
  # 3. 将 Buffer 转为 sf 对象（exact_extract 的要求）
  buffer_sf <- st_as_sf(buffer_vec)
  
  # 4. 执行分区统计提取
  message("正在提取缓冲区统计信息...")
  
  # 海拔和坡度提取平均值 (mean)
  stats_continuous <- exact_extract(terrain_stack[[c("elev", "slope")]], 
                                    buffer_sf, 
                                    fun = "mean", 
                                    progress = FALSE)
  
  # 坡向提取众数 (majority) - 专门定义的众数逻辑
  # 避免 1度和359度平均成180度的问题
  aspect_mode <- exact_extract(terrain_stack[["aspect"]], 
                               buffer_sf, 
                               function(values, coverage_fraction) {
                                 v <- values[values >= 0 & !is.na(values)]
                                 if(length(v) == 0) return(NA)
                                 ux <- unique(v)
                                 ux[which.max(tabulate(match(v, ux)))]
                               }, 
                               progress = FALSE)
  
  # 5. 组装结果
  result <- buffer_sf %>%
    mutate(
      mean_elev = stats_continuous$mean.elev,
      mean_slope = stats_continuous$mean.slope,
      mode_aspect = aspect_mode
    )
  
  message("提取完成！")
  return(result)
}









calculate_lulc_loss_1_34 <- function(lulc_2002, lulc_2018, buffer_vec) {
  
  # 1. 确保栅格对齐 (如果 2018 和 2002 不对齐，先进行重采样)
  if (!compareGeom(lulc_2002, lulc_2018, stopOnError = FALSE)) {
    message("栅格不对齐，正在重采样 lulc_2018...")
    lulc_2018 <- resample(lulc_2018, lulc_2002, method = "near")
  }
  
  # 2. 识别变化像素
  # 逻辑：如果 2002 为 1 且 2018 为 3,4,5，则标记为 1，否则为 0
  # 使用 ifel 函数处理大规模栅格非常高效
  change_layer <- ifel(lulc_2002 == 1 & (lulc_2018 == 3), 1, 0)
  #change_layer <- ifel((lulc_2002 == 1 | lulc_2002 == 2) & (lulc_2018 == 3 | lulc_2018 == 4 | lulc_2018 == 5), 1, 0)
  
  # 3. 计算单个像素的面积 (单位：平方米)
  # res(change_layer)[1] 获取横向分辨率，[2] 获取纵向
  pixel_area <- res(change_layer)[1] * res(change_layer)[2]
  
  # 4. 在 Buffer 内统计变化像素的数量
  # sum 会计算 buffer 范围内所有值为 1 的像素的总和
  buffer_sf <- st_as_sf(buffer_vec)
  change_counts <- exact_extract(change_layer, buffer_sf, 'sum', progress = TRUE)
  
  # 5. 换算为面积（单位：平方米，如果需要平方公里则再除以 1,000,000）
  loss_area_m2 <- change_counts * pixel_area
  
  return(loss_area_m2)
}







calculate_lulc_loss_1_5 <- function(lulc_2002, lulc_2018, buffer_vec) {
  
  # 1. 确保栅格对齐 (如果 2018 和 2002 不对齐，先进行重采样)
  if (!compareGeom(lulc_2002, lulc_2018, stopOnError = FALSE)) {
    message("栅格不对齐，正在重采样 lulc_2018...")
    lulc_2018 <- resample(lulc_2018, lulc_2002, method = "near")
  }
  
  # 2. 识别变化像素
  # 逻辑：如果 2002 为 1 且 2018 为 3,4,5，则标记为 1，否则为 0
  # 使用 ifel 函数处理大规模栅格非常高效
  change_layer <- ifel(lulc_2002 == 1 & (lulc_2018 == 5), 1, 0)
  #change_layer <- ifel((lulc_2002 == 1 | lulc_2002 == 2) & (lulc_2018 == 3 | lulc_2018 == 4 | lulc_2018 == 5), 1, 0)
  
  # 3. 计算单个像素的面积 (单位：平方米)
  # res(change_layer)[1] 获取横向分辨率，[2] 获取纵向
  pixel_area <- res(change_layer)[1] * res(change_layer)[2]
  
  # 4. 在 Buffer 内统计变化像素的数量
  # sum 会计算 buffer 范围内所有值为 1 的像素的总和
  buffer_sf <- st_as_sf(buffer_vec)
  change_counts <- exact_extract(change_layer, buffer_sf, 'sum', progress = TRUE)
  
  # 5. 换算为面积（单位：平方米，如果需要平方公里则再除以 1,000,000）
  loss_area_m2 <- change_counts * pixel_area
  
  return(loss_area_m2)
}





calculate_lulc_loss_12_5 <- function(lulc_2002, lulc_2018, buffer_vec) {
  
  # 1. 确保栅格对齐 (如果 2018 和 2002 不对齐，先进行重采样)
  if (!compareGeom(lulc_2002, lulc_2018, stopOnError = FALSE)) {
    message("栅格不对齐，正在重采样 lulc_2018...")
    lulc_2018 <- resample(lulc_2018, lulc_2002, method = "near")
  }
  
  # 2. 识别变化像素
  # 逻辑：如果 2002 为 1 且 2018 为 3,4,5，则标记为 1，否则为 0
  # 使用 ifel 函数处理大规模栅格非常高效
  #change_layer <- ifel(lulc_2002 == 1 & (lulc_2018 == 3 | lulc_2018 == 4 | lulc_2018 == 5), 1, 0)
  change_layer <- ifel((lulc_2002 == 2) & (lulc_2018 == 5), 1, 0)
  
  # 3. 计算单个像素的面积 (单位：平方米)
  # res(change_layer)[1] 获取横向分辨率，[2] 获取纵向
  pixel_area <- res(change_layer)[1] * res(change_layer)[2]
  
  # 4. 在 Buffer 内统计变化像素的数量
  # sum 会计算 buffer 范围内所有值为 1 的像素的总和
  buffer_sf <- st_as_sf(buffer_vec)
  change_counts <- exact_extract(change_layer, buffer_sf, 'sum', progress = TRUE)
  
  # 5. 换算为面积（单位：平方米，如果需要平方公里则再除以 1,000,000）
  loss_area_m2 <- change_counts * pixel_area
  
  return(loss_area_m2)
}




calculate_lulc_loss_12_34 <- function(lulc_2002, lulc_2018, buffer_vec) {
  
  # 1. 确保栅格对齐 (如果 2018 和 2002 不对齐，先进行重采样)
  if (!compareGeom(lulc_2002, lulc_2018, stopOnError = FALSE)) {
    message("栅格不对齐，正在重采样 lulc_2018...")
    lulc_2018 <- resample(lulc_2018, lulc_2002, method = "near")
  }
  
  # 2. 识别变化像素
  # 逻辑：如果 2002 为 1 且 2018 为 3,4,5，则标记为 1，否则为 0
  # 使用 ifel 函数处理大规模栅格非常高效
  #change_layer <- ifel(lulc_2002 == 1 & (lulc_2018 == 3 | lulc_2018 == 4 | lulc_2018 == 5), 1, 0)
  change_layer <- ifel((lulc_2002 == 2) & (lulc_2018 == 3), 1, 0)
  
  # 3. 计算单个像素的面积 (单位：平方米)
  # res(change_layer)[1] 获取横向分辨率，[2] 获取纵向
  pixel_area <- res(change_layer)[1] * res(change_layer)[2]
  
  # 4. 在 Buffer 内统计变化像素的数量
  # sum 会计算 buffer 范围内所有值为 1 的像素的总和
  buffer_sf <- st_as_sf(buffer_vec)
  change_counts <- exact_extract(change_layer, buffer_sf, 'sum', progress = TRUE)
  
  # 5. 换算为面积（单位：平方米，如果需要平方公里则再除以 1,000,000）
  loss_area_m2 <- change_counts * pixel_area
  
  return(loss_area_m2)
}








calculate_lulc_area <- function(lulc_raster, buffer_vec, class_vals) {
  
  # 1. 把目标地类变成二值栅格：目标类 = 1，其余 = 0
  class_mask <- Reduce(`|`, lapply(class_vals, function(v) lulc_raster == v))
  area_layer <- ifel(class_mask, 1, 0)
  
  # 2. 计算单个像元面积（平方米）
  pixel_area <- res(area_layer)[1] * res(area_layer)[2]
  
  # 3. 转成 sf，供 exact_extract 使用
  buffer_sf <- st_as_sf(buffer_vec)
  
  # 4. 统计每个 buffer 内目标像元数量
  class_counts <- exactextractr::exact_extract(
    area_layer,
    buffer_sf,
    "sum",
    progress = FALSE
  )
  
  # 5. 像元数量 × 单像元面积 = 目标地类面积（平方米）
  area_m2 <- class_counts * pixel_area
  
  return(area_m2)
}











#' 样点分层筛选函数
#' @param study_area sf对象，撒点的空间范围（如村庄边界）
#' @param lulc_raster SpatRaster对象，用于筛选初始地类的栅格
#' @param target_lulc 向量，允许保留的 LULC 类别代码（如 c(1, 2)）
#' @param exclude_polygons_list 列表，包含所有需要剔除的 sf 对象（如 list(PA, Buffer)）
#' @param initial_size 整数，初始撒点的数量
#' @param final_size 整数，最终期望保留的数量
#' @param seed 整数，随机种子

select_stratified_points <- function(study_area, 
                                     lulc_raster, 
                                     target_lulc = c(1), 
                                     exclude_polygons_list = list(), 
                                     initial_size = 30000, 
                                     final_size = 10000) {
  
  # 1. 在研究区域内随机撒点
  message("正在生成初始随机点...")
  pts <- st_sample(study_area, size = initial_size) %>% st_as_sf()
  
  # 2. 筛选特定的 LULC 类型
  message("正在根据 LULC 类别进行第一轮筛选...")
  # 提取栅格值（ID列通常是第一列，值是第二列）
  extracted_vals <- terra::extract(lulc_raster, vect(pts))
  pts$lulc_val <- extracted_vals[, 2]
  
  pts <- pts[pts$lulc_val %in% target_lulc, ]
  
  # 3. 循环剔除不需要的区域 (PA, Buffers 等)
  if (length(exclude_polygons_list) > 0) {
    message("正在根据排除区域进行筛选...")
    for (i in seq_along(exclude_polygons_list)) {
      # 找出落在排除区内的点
      is_inside <- st_intersects(pts, exclude_polygons_list[[i]], sparse = FALSE) %>% 
        apply(1, any)
      # 只保留不在排除区内的点
      pts <- pts[!is_inside, ]
    }
  }
  
  # 4. 最后抽样
  current_count <- nrow(pts)
  if (current_count >= final_size) {
    message("筛选完成，正在提取最终样点...")
    pts <- pts[sample(1:current_count, final_size), ]
  } else {
    warning("合格样点不足 ", final_size, " 个，当前仅剩下 ", current_count, " 个。")
  }
  
  return(pts)
}




library(sf)

#' 计算最近邻距离并返回等长向量
#' @param pts_A 原始观测点 (sf 对象)
#' @param pts_B 参照点 (sf 对象)
#' @param range_poly 限制范围的多边形 (sf 对象)
#' @return 返回一个数值向量，长度与 pts_A 相同

get_dist_vector <- function(pts_A, pts_B, range_poly) {
  
  # 1. 统一坐标系（以多边形为准）
  target_crs <- st_crs(range_poly)
  pts_A_trans <- st_transform(pts_A, target_crs)
  pts_B_trans <- st_transform(pts_B, target_crs)
  
  # 2. 识别哪些点在多边形内 (返回逻辑向量)
  # sparse = FALSE 得到 TRUE/FALSE 向量
  is_inside <- lengths(st_intersects(pts_A_trans, range_poly)) > 0
  
  # 3. 初始化结果向量（默认填充 NA）
  dist_vec <- rep(NA_real_, nrow(pts_A))
  
  # 4. 仅对范围内的点计算最近邻
  if (any(is_inside)) {
    # 提取范围内的点
    pts_A_in <- pts_A_trans[is_inside, ]
    
    # 寻找最近邻索引
    nn_idx <- st_nearest_feature(pts_A_in, pts_B_trans)
    
    # 计算距离并存入对应位置
    dist_vec[is_inside] <- as.numeric(
      st_distance(pts_A_in, pts_B_trans[nn_idx, ], by_element = TRUE)
    )
  }
  
  return(dist_vec)
}