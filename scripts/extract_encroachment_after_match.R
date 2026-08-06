


#pseudo在1/2/3/4/5，initial forest land不为0

matched_forest_dai_hani$loss_1_34_14 <- calculate_lulc_loss_1_34(lulc_2002, lulc_2014, matched_forest_dai_hani)
matched_forest_dai_hani$loss_1_34_intensity_14 <- matched_forest_dai_hani$loss_1_34_14/matched_forest_dai_hani$forest_area_02

matched_forest_dai_hani$loss_1_345_14 <- calculate_lulc_loss_1_5(lulc_2002, lulc_2014, matched_forest_dai_hani)
matched_forest_dai_hani$loss_1_345_intensity_14 <- matched_forest_dai_hani$loss_1_345_14/matched_forest_dai_hani$forest_area_02



#pseudo在1/2/3/4/5，initial woody land不为0

# st_write(final_buffer_500, "final_buffer_500.gpkg", delete_dsn = TRUE)


matched_woody_dai_hani$loss_12_34_14 <- calculate_lulc_loss_12_34(lulc_2002, lulc_2014, matched_woody_dai_hani)
matched_woody_dai_hani$loss_12_34_intensity_14 <- matched_woody_dai_hani$loss_12_34_14/matched_woody_dai_hani$woody_area_02

matched_woody_dai_hani$loss_12_345_14 <- calculate_lulc_loss_12_5(lulc_2002, lulc_2014, matched_woody_dai_hani)
matched_woody_dai_hani$loss_12_345_intensity_14 <- matched_woody_dai_hani$loss_12_345_14/matched_woody_dai_hani$woody_area_02

if(FALSE){

#pseudo在1/2/3/4/5，initial forest land不为0

final_buffer_500_12345_forest$loss_1_34_14 <- calculate_lulc_loss_1_34(lulc_2002, lulc_2014, final_buffer_500_12345_forest)
final_buffer_500_12345_forest$loss_1_34_intensity_14 <- final_buffer_500_12345_forest$loss_1_34_14/final_buffer_500_12345_forest$forest_area_02

final_buffer_500_12345_forest$loss_1_345_14 <- calculate_lulc_loss_1_5(lulc_2002, lulc_2014, final_buffer_500_12345_forest)
final_buffer_500_12345_forest$loss_1_345_intensity_14 <- final_buffer_500_12345_forest$loss_1_345_14/final_buffer_500_12345_forest$forest_area_02



#pseudo在1/2/3/4/5，initial woody land不为0

# st_write(final_buffer_500, "final_buffer_500.gpkg", delete_dsn = TRUE)


final_buffer_500_12345_woody$loss_12_34_14 <- calculate_lulc_loss_12_34(lulc_2002, lulc_2014, final_buffer_500_12345_woody)
final_buffer_500_12345_woody$loss_12_34_intensity_14 <- final_buffer_500_12345_woody$loss_12_34_14/final_buffer_500_12345_woody$woody_area_02

final_buffer_500_12345_woody$loss_12_345_14 <- calculate_lulc_loss_12_5(lulc_2002, lulc_2014, final_buffer_500_12345_woody)
final_buffer_500_12345_woody$loss_12_345_intensity_14 <- final_buffer_500_12345_woody$loss_12_345_14/final_buffer_500_12345_woody$woody_area_02

}

