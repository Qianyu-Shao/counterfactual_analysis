##########################
# pseudo from all land, buffer has forest, forest-cashcrop

# 2002-2014

matched_forest_dai_hani <- matched_forest_dai_hani  %>%
  mutate(log_income = log(net_income_adjusted),
         log_popu = log(popu_density),
         log_agri = log(agri_land_pp))

# 先跑一个普通 OLS
ols_fit <- lm(loss_1_34_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
              data = matched_forest_dai_hani )

summary(ols_fit)
vif(ols_fit)

centroids <- st_centroid(matched_forest_dai_hani [names(residuals(ols_fit)), ])
coords_pure <- as.matrix(st_coordinates(centroids)[, 1:2])

knn_obj <- knearneigh(coords_pure, k = 5)

# 5. 将 kNN 转换为邻接列表 (nb)
nb_knn <- knn2nb(knn_obj)

# 6. 将邻接列表转换为行标准化权重矩阵 (listw)
# style = "W" 代表行标准化，这是回归分析的标准要求
W_knn <- nb2listw(nb_knn, style = "W")


rs_results <- lm.RStests(ols_fit, W_knn, test = "all")
summary(rs_results)





# 拟合 SLM 模型 (SAR)
slm_fit_1_3 <- spatialreg::lagsarlm(
  loss_1_34_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * hani, 
  data = matched_forest_dai_hani , 
  listw = W_knn
)


summary(slm_fit_1_3)



##########################

# pseudo from all land, buffer has forest, forest-farmland


# 先跑一个普通 OLS
ols_fit <- lm(loss_1_345_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
              data = matched_forest_dai_hani )


rs_results <- lm.RStests(ols_fit, W_knn, test = "all")
summary(rs_results)



# 拟合 SLM 模型 (SAR)
slm_fit_1_5 <- spatialreg::lagsarlm(
  loss_1_345_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * hani, 
  data = matched_forest_dai_hani , 
  listw = W_knn
)


summary(slm_fit_1_5)



##########################

# pseudo from all land, buffer has woody, woody to cashcrop

matched_woody_dai_hani  <- matched_woody_dai_hani  %>%
  mutate(log_income = log(net_income_adjusted),
         log_popu = log(popu_density),
         log_agri = log(agri_land_pp))


# 先跑一个普通 OLS
ols_fit <- lm(loss_12_34_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
              data = matched_woody_dai_hani )

vif(ols_fit)

centroids <- st_centroid(matched_woody_dai_hani [names(residuals(ols_fit)), ])
coords_pure <- as.matrix(st_coordinates(centroids)[, 1:2])

knn_obj <- knearneigh(coords_pure, k = 5)

# 5. 将 kNN 转换为邻接列表 (nb)
nb_knn <- knn2nb(knn_obj)

# 6. 将邻接列表转换为行标准化权重矩阵 (listw)
# style = "W" 代表行标准化，这是回归分析的标准要求
W_knn <- nb2listw(nb_knn, style = "W")


rs_results <- lm.RStests(ols_fit, W_knn, test = "all")
summary(rs_results)



# 拟合 SLM 模型 (SAR)
slm_fit_12_3 <- spatialreg::lagsarlm(
  loss_12_34_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
  data = matched_woody_dai_hani , 
  listw = W_knn
)


summary(slm_fit_12_3)


# st_write(matched_woody_dai_hani, "matched_woody_dai_hani.gpkg", delete_dsn = TRUE)

##########################


## pseudo from all land, buffer has woody, woody-farmland


# 先跑一个普通 OLS
ols_fit <- lm(loss_12_345_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
              data = matched_woody_dai_hani )



rs_results <- lm.RStests(ols_fit, W_knn, test = "all")
summary(rs_results)




# 拟合 SLM 模型 (SAR)
slm_fit_12_5 <- spatialreg::lagsarlm(
  loss_12_345_intensity_14 ~ log_income + log_agri + log_popu + is_sns + is_sns * dai, 
  data = matched_woody_dai_hani , 
  listw = W_knn
)


summary(slm_fit_12_5)
























