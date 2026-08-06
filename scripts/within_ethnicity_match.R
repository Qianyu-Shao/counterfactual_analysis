library(dplyr)
library(MatchIt)
library(sf)

# 如果是 sf 对象，MatchIt 一般更稳的是先 drop geometry
forest_df <- final_buffer_500_12345_forest %>%
  mutate(
    ethnicity = case_when(
      dai == 1 ~ "Dai",
      hani == 1 ~ "Hani",
      TRUE ~ "Bulang"
    )
  ) %>%
  filter(ethnicity %in% c("Dai", "Hani", "Bulang"))

woody_df <- final_buffer_500_12345_woody %>%
  mutate(
    ethnicity = case_when(
      dai == 1 ~ "Dai",
      hani == 1 ~ "Hani",
      TRUE ~ "Bulang"
    )
  ) %>%
  filter(ethnicity %in% c("Dai", "Hani", "Bulang"))



run_match_by_ethnicity_forest <- function(data, eth_name, caliper = 0.05) {
  subdat <- data %>%
    filter(ethnicity == eth_name)
  
  m <- matchit(
    is_sns ~ mean_elev + mean_slope + dist_to_road + sin(mode_aspect * pi / 180) + cos(mode_aspect * pi / 180) + dist_to_village,
#    is_sns ~ mean_elev + mean_slope + sin(mode_aspect * pi / 180) + cos(mode_aspect * pi / 180),
#    is_sns ~ mean_elev + mean_slope + dist_to_road + mode_aspect + forest_area_02 + dist_to_village,
    data = subdat,
    method = "nearest",
    distance = "glm",
    caliper = caliper,
    replace = FALSE
  )
  
  matched <- match.data(m) %>%
    mutate(ethnicity = eth_name)
  
  list(
    match_obj = m,
    matched_data = matched
  )
}


run_match_by_ethnicity_woody <- function(data, eth_name, caliper = 0.05) {
  subdat <- data %>%
    filter(ethnicity == eth_name)
  
  m <- matchit(
#    is_sns ~ mean_elev + mean_slope + dist_to_road + mode_aspect + woody_area_02 + dist_to_village,
   is_sns ~ mean_elev + mean_slope + dist_to_road + cos(mode_aspect * pi / 180) + sin(mode_aspect * pi / 180) + dist_to_village,
#    is_sns ~ mean_elev + mean_slope + cos(mode_aspect * pi / 180) + sin(mode_aspect * pi / 180),
    data = subdat,
    method = "nearest",
    distance = "glm",
    caliper = caliper,
    replace = FALSE
  )
  
  matched <- match.data(m) %>%
    mutate(ethnicity = eth_name)
  
  list(
    match_obj = m,
    matched_data = matched
  )
}


forest_match_dai <- run_match_by_ethnicity_forest(forest_df, "Dai", caliper = 0.03)
forest_match_hani <- run_match_by_ethnicity_forest(forest_df, "Hani", caliper = 0.05)
forest_match_bulang <- run_match_by_ethnicity_forest(forest_df, "Bulang", caliper = 0.05)



woody_match_dai <- run_match_by_ethnicity_woody(woody_df, "Dai", caliper = 0.03)
woody_match_hani <- run_match_by_ethnicity_woody(woody_df, "Hani", caliper = 0.05)
woody_match_bulang <- run_match_by_ethnicity_woody(woody_df, "Bulang", caliper = 0.05)




matched_forest_dai <- forest_match_dai$matched_data
matched_forest_hani <- forest_match_hani$matched_data
matched_forest_bulang <- forest_match_bulang$matched_data



matched_forest_all_ethnicity <- bind_rows(
  matched_forest_dai,
  matched_forest_hani,
  matched_forest_bulang
)

matched_forest_dai_hani <- bind_rows(
  matched_forest_dai,
  matched_forest_hani
)




matched_woody_dai <- woody_match_dai$matched_data
matched_woody_hani <- woody_match_hani$matched_data
matched_woody_bulang <- woody_match_bulang$matched_data

matched_woody_all_ethnicity <- bind_rows(
  matched_woody_dai,
  matched_woody_hani,
  matched_woody_bulang
)

matched_woody_dai_hani <- bind_rows(
  matched_woody_dai,
  matched_woody_hani
)


table(matched_forest_all_ethnicity$ethnicity, matched_forest_all_ethnicity$is_sns)
table(matched_woody_all_ethnicity$ethnicity, matched_woody_all_ethnicity$is_sns)



summary(forest_match_dai$match_obj)
summary(forest_match_hani$match_obj)
summary(forest_match_bulang$match_obj)

summary(woody_match_dai$match_obj)
summary(woody_match_hani$match_obj)
summary(woody_match_bulang$match_obj)


# -------------------------------------------------------
# Export concise matching diagnostics for paper reporting
# -------------------------------------------------------

extract_match_counts <- function(data, matched_data, sample_name, eth_name) {
  before_counts <- data %>%
    filter(ethnicity == eth_name) %>%
    count(is_sns, name = "n") %>%
    mutate(stage = "Before")

  after_counts <- matched_data %>%
    count(is_sns, name = "n") %>%
    mutate(stage = "After")

  bind_rows(before_counts, after_counts) %>%
    mutate(
      sample = sample_name,
      ethnicity = eth_name,
      site_type = ifelse(is_sns == 1, "True SNS", "Pseudo-SNS")
    ) %>%
    select(sample, ethnicity, stage, site_type, n)
}

extract_balance_table <- function(match_obj, sample_name, eth_name) {
  smry <- summary(match_obj, standardize = TRUE)

  all_tab <- as.data.frame(smry$sum.all)
  all_tab$covariate <- rownames(all_tab)
  rownames(all_tab) <- NULL

  matched_tab <- as.data.frame(smry$sum.matched)
  matched_tab$covariate <- rownames(matched_tab)
  rownames(matched_tab) <- NULL

  all_smd_col <- grep("Std\\.", names(all_tab), value = TRUE)[1]
  matched_smd_col <- grep("Std\\.", names(matched_tab), value = TRUE)[1]

  dplyr::full_join(
    all_tab %>%
      transmute(
        covariate,
        abs_smd_before = abs(.data[[all_smd_col]])
      ),
    matched_tab %>%
      transmute(
        covariate,
        abs_smd_after = abs(.data[[matched_smd_col]])
      ),
    by = "covariate"
  ) %>%
    mutate(
      sample = sample_name,
      ethnicity = eth_name,
      improved = abs_smd_after < abs_smd_before
    ) %>%
    select(sample, ethnicity, covariate, abs_smd_before, abs_smd_after, improved) %>%
    mutate(
      abs_smd_before = round(abs_smd_before, 3),
      abs_smd_after = round(abs_smd_after, 3)
    )
}

matching_counts_table <- bind_rows(
  extract_match_counts(forest_df, matched_forest_dai, "Forest", "Dai"),
  extract_match_counts(forest_df, matched_forest_hani, "Forest", "Hani"),
  extract_match_counts(forest_df, matched_forest_bulang, "Forest", "Bulang"),
  extract_match_counts(woody_df, matched_woody_dai, "Woody", "Dai"),
  extract_match_counts(woody_df, matched_woody_hani, "Woody", "Hani"),
  extract_match_counts(woody_df, matched_woody_bulang, "Woody", "Bulang")
)

matching_balance_table <- bind_rows(
  extract_balance_table(forest_match_dai$match_obj, "Forest", "Dai"),
  extract_balance_table(forest_match_hani$match_obj, "Forest", "Hani"),
  extract_balance_table(forest_match_bulang$match_obj, "Forest", "Bulang"),
  extract_balance_table(woody_match_dai$match_obj, "Woody", "Dai"),
  extract_balance_table(woody_match_hani$match_obj, "Woody", "Hani"),
  extract_balance_table(woody_match_bulang$match_obj, "Woody", "Bulang")
)

dir.create("./plot", showWarnings = FALSE, recursive = TRUE)

write.csv(
  matching_counts_table,
  "./plot/within_ethnicity_matching_counts.csv",
  row.names = FALSE
)

write.csv(
  matching_balance_table,
  "./plot/within_ethnicity_matching_balance.csv",
  row.names = FALSE
)

matching_maintext_table <- matching_counts_table %>%
  filter(ethnicity %in% c("Dai", "Hani")) %>%
  tidyr::pivot_wider(
    names_from = c(stage, site_type),
    values_from = n
  ) %>%
  left_join(
    matching_balance_table %>%
      filter(ethnicity %in% c("Dai", "Hani")) %>%
      group_by(sample, ethnicity) %>%
      summarise(
        mean_abs_smd_before = mean(abs_smd_before, na.rm = TRUE),
        mean_abs_smd_after = mean(abs_smd_after, na.rm = TRUE),
        max_abs_smd_before = max(abs_smd_before, na.rm = TRUE),
        max_abs_smd_after = max(abs_smd_after, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("sample", "ethnicity")
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  ) %>%
  arrange(sample, ethnicity)

write.csv(
  matching_maintext_table,
  "./plot/within_ethnicity_matching_maintext.csv",
  row.names = FALSE
)

message("Saved concise matching diagnostics:")
message("./plot/within_ethnicity_matching_counts.csv")
message("./plot/within_ethnicity_matching_balance.csv")
message("./plot/within_ethnicity_matching_maintext.csv")


