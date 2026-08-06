# ============================================================
# Dai vs Hani: Geography vs Institution
# Purpose:
# 1. Visualize current geographic differences between Dai and Hani true sacred sites
# 2. Test whether Dai-Hani differences in SNS buffering are mainly due to
#    geography or remain after geography is made comparable
# 3. Provide a stepwise analysis workflow with starter code
# ============================================================

# ----------------------------
# 0. Load the current pipeline
# ----------------------------
# Assumption:
# - function.R defines helper functions
# - dai_with_temple.R creates the main pseudo/true datasets with geography
# - within_ethnicity_match.R creates matched_forest_dai_hani / matched_woody_dai_hani
# - extract_encroachment_after_match.R adds encroachment variables
#
# Run this file in a fresh session if possible.

source("./scripts/function.R")
source("./scripts/bulang_clean.R")
source("./scripts/within_ethnicity_match.R")
source("./scripts/extract_encroachment_after_match.R")

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(spdep)
  library(spatialreg)
  library(rstatix)
  library(MatchIt)
})

# -------------------------------------------------------
# 1. Start with true sacred sites only: geography profile
# -------------------------------------------------------
# Goal:
# Show whether Dai and Hani true sacred sites are already located
# in very different terrain / accessibility environments.

make_aspect_sector <- function(x) {
  case_when(
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

true_dh_forest <- matched_forest_dai_hani %>%
  filter(is_sns == 1, dai == 1 | hani == 1) %>%
  mutate(
    ethnicity = ifelse(dai == 1, "Dai", "Hani"),
    baseline = "Forest",
    aspect_sector = make_aspect_sector(mode_aspect)
  )

true_dh_woody <- matched_woody_dai_hani %>%
  filter(is_sns == 1, dai == 1 | hani == 1) %>%
  mutate(
    ethnicity = ifelse(dai == 1, "Dai", "Hani"),
    baseline = "Woody",
    aspect_sector = make_aspect_sector(mode_aspect)
  )

true_dh_all <- bind_rows(true_dh_forest, true_dh_woody) %>%
  mutate(
    ethnicity = factor(ethnicity, levels = c("Dai", "Hani")),
    baseline = factor(baseline, levels = c("Forest", "Woody"))
  )

plot_density_with_points <- function(df, baseline_name, file_name) {
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

  point_df <- density_df %>%
    group_by(ethnicity, variable) %>%
    group_modify(~ {
      d <- density(.x$value, na.rm = TRUE)
      .x$density_y <- approx(
        x = d$x,
        y = d$y,
        xout = .x$value,
        rule = 2
      )$y
      .x
    }) %>%
    ungroup()

  p <- ggplot(
    density_df,
    aes(x = value, color = ethnicity, fill = ethnicity)
  ) +
    geom_density(alpha = 0.12, adjust = 1.1, linewidth = 0.9) +
    geom_point(
      data = point_df,
      aes(y = density_y, shape = ethnicity),
      alpha = 0.65,
      size = 1.6,
      stroke = 0,
      inherit.aes = TRUE
    ) +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    scale_color_manual(values = c("Dai" = "#D95F02", "Hani" = "#1B9E77")) +
    scale_fill_manual(values = c("Dai" = "#D95F02", "Hani" = "#1B9E77")) +
    scale_shape_manual(values = c("Dai" = 16, "Hani" = 17)) +
    labs(
      title = paste("Geographic Distributions of Dai and Hani True Sacred Sites:", baseline_name),
      x = NULL,
      y = "Density",
      color = NULL,
      fill = NULL,
      shape = NULL
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

p_density_forest <- plot_density_with_points(
  true_dh_forest,
  "Forest",
  "./plot/dai_hani_true_sites_density_forest.png"
)

p_density_woody <- plot_density_with_points(
  true_dh_woody,
  "Woody",
  "./plot/dai_hani_true_sites_density_woody.png"
)

# 1B. Aspect rose plots
aspect_df <- true_dh_all %>%
  st_drop_geometry() %>%
  count(baseline, ethnicity, aspect_sector) %>%
  filter(!is.na(aspect_sector)) %>%
  mutate(
    aspect_sector = factor(
      aspect_sector,
      levels = c("North", "Northeast", "East", "Southeast",
                 "South", "Southwest", "West", "Northwest")
    )
  )

plot_aspect_by_baseline <- function(df, baseline_name, file_name) {
  p <- ggplot(df, aes(x = aspect_sector, y = n, fill = ethnicity)) +
    geom_col(alpha = 0.95) +
    coord_polar() +
    facet_wrap(~ ethnicity, nrow = 1) +
    scale_fill_manual(values = c("Dai" = "#D95F02", "Hani" = "#1B9E77")) +
    labs(
      title = paste("Aspect Distribution of Dai and Hani True Sacred Sites:", baseline_name),
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

p_aspect_forest <- plot_aspect_by_baseline(
  aspect_df %>% filter(baseline == "Forest"),
  "Forest",
  "./plot/dai_hani_true_sites_aspect_forest.png"
)

p_aspect_woody <- plot_aspect_by_baseline(
  aspect_df %>% filter(baseline == "Woody"),
  "Woody",
  "./plot/dai_hani_true_sites_aspect_woody.png"
)

# 1C. Quick numeric comparison on true sacred sites
true_geo_compare <- true_dh_all %>%
  st_drop_geometry() %>%
  select(ethnicity, baseline, mean_elev, mean_slope, dist_to_road, dist_to_village) %>%
  pivot_longer(
    cols = c(mean_elev, mean_slope, dist_to_road, dist_to_village),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(baseline, variable) %>%
  wilcox_test(value ~ ethnicity) %>%
  ungroup()

write.csv(
  true_geo_compare,
  "./plot/dai_hani_true_sites_geography_compare.csv",
  row.names = FALSE
)

# ---------------------------------------------------
# 2. Build a Dai-Hani common-support sample (overlap)
# ---------------------------------------------------
# Goal:
# Restrict comparison to locations where Dai and Hani observations
# occupy comparable geography. If ethnic SNS differences disappear
# here, geography likely explains much of the gap.

prepare_overlap_sample <- function(data) {
  data %>%
    filter(dai == 1 | hani == 1) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

forest_dh <- prepare_overlap_sample(matched_forest_dai_hani)
woody_dh  <- prepare_overlap_sample(matched_woody_dai_hani)

build_ethnicity_overlap <- function(df) {
  base_df <- df %>%
    st_drop_geometry() %>%
    mutate(
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )

  ps_fit <- glm(
    as.formula(
      paste(
        "dai_eth ~ mean_elev + mean_slope + aspect_cos + aspect_sin +",
        "dist_to_road + dist_to_village"
      )
    ),
    family = binomial(),
    data = base_df
  )

  base_df$ps_eth <- predict(ps_fit, type = "response")

  dai_range  <- range(base_df$ps_eth[base_df$dai_eth == 1], na.rm = TRUE)
  hani_range <- range(base_df$ps_eth[base_df$dai_eth == 0], na.rm = TRUE)

  lower <- max(dai_range[1], hani_range[1])
  upper <- min(dai_range[2], hani_range[2])

  overlap_ids <- base_df %>%
    filter(ps_eth >= lower, ps_eth <= upper) %>%
    pull(ID)

  df_overlap <- df %>%
    filter(ID %in% overlap_ids) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )

  list(
    model = ps_fit,
    overlap_df = df_overlap,
    overlap_bounds = c(lower = lower, upper = upper),
    ps_data = base_df
  )
}

forest_overlap <- build_ethnicity_overlap(forest_dh)
woody_overlap  <- build_ethnicity_overlap(woody_dh)

# Optional: visualize PS overlap
plot_ps_overlap <- function(ps_data, title, file_name) {
  p <- ggplot(ps_data, aes(x = ps_eth, fill = factor(dai_eth), color = factor(dai_eth))) +
    geom_density(alpha = 0.18) +
    scale_fill_manual(values = c("0" = "#1B9E77", "1" = "#D95F02"),
                      labels = c("Hani", "Dai")) +
    scale_color_manual(values = c("0" = "#1B9E77", "1" = "#D95F02"),
                       labels = c("Hani", "Dai")) +
    labs(title = title, x = "Ethnicity propensity score", y = "Density", fill = NULL, color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(file_name, p, width = 8, height = 5, dpi = 300)
}

plot_ps_overlap(
  forest_overlap$ps_data,
  "Dai-Hani Common Support: Forest Sample",
  "./plot/dai_hani_overlap_forest_ps.png"
)

plot_ps_overlap(
  woody_overlap$ps_data,
  "Dai-Hani Common Support: Woody Sample",
  "./plot/dai_hani_overlap_woody_ps.png"
)





check_overlap_counts <- function(full_df, overlap_df, sample_name) {
  before_tab <- full_df %>%
    st_drop_geometry() %>%
    mutate(
      ethnicity = ifelse(dai == 1, "Dai", "Hani"),
      stage = "Before overlap"
    ) %>%
    count(stage, ethnicity, name = "n")
  
  after_tab <- overlap_df %>%
    st_drop_geometry() %>%
    mutate(
      ethnicity = ifelse(dai == 1, "Dai", "Hani"),
      stage = "After overlap"
    ) %>%
    count(stage, ethnicity, name = "n")
  
  bind_rows(before_tab, after_tab) %>%
    mutate(sample = sample_name) %>%
    select(sample, stage, ethnicity, n)
}

forest_overlap_counts <- check_overlap_counts(
  forest_dh,
  forest_overlap$overlap_df,
  "Forest"
)

woody_overlap_counts <- check_overlap_counts(
  woody_dh,
  woody_overlap$overlap_df,
  "Woody"
)

overlap_counts_all <- bind_rows(forest_overlap_counts, woody_overlap_counts)

print(overlap_counts_all)





calc_smd_two_groups <- function(data, var, group_var = "dai_eth") {
  x1 <- data[[var]][data[[group_var]] == 1]
  x0 <- data[[var]][data[[group_var]] == 0]
  
  m1 <- mean(x1, na.rm = TRUE)
  m0 <- mean(x0, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE)
  s0 <- sd(x0, na.rm = TRUE)
  
  pooled_sd <- sqrt((s1^2 + s0^2) / 2)
  smd <- ifelse(is.na(pooled_sd) || pooled_sd == 0, NA, (m1 - m0) / pooled_sd)
  
  data.frame(
    variable = var,
    mean_dai = m1,
    mean_hani = m0,
    smd = smd
  )
}

check_balance_before_after_overlap <- function(full_df, overlap_df, vars, sample_name) {
  full_clean <- full_df %>%
    st_drop_geometry() %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )
  
  overlap_clean <- overlap_df %>%
    st_drop_geometry() %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )
  
  before_res <- bind_rows(lapply(vars, function(v) {
    calc_smd_two_groups(full_clean, v) %>%
      mutate(stage = "Before overlap")
  }))
  
  after_res <- bind_rows(lapply(vars, function(v) {
    calc_smd_two_groups(overlap_clean, v) %>%
      mutate(stage = "After overlap")
  }))
  
  bind_rows(before_res, after_res) %>%
    mutate(sample = sample_name) %>%
    select(sample, stage, everything())
}


forest_balance_overlap <- check_balance_before_after_overlap(
  forest_dh,
  forest_overlap$overlap_df,
  vars = c("mean_elev", "mean_slope", "aspect_cos", "aspect_sin", "dist_to_road", "dist_to_village"),
  sample_name = "Forest"
)

woody_balance_overlap <- check_balance_before_after_overlap(
  woody_dh,
  woody_overlap$overlap_df,
  vars = c("mean_elev", "mean_slope", "aspect_cos", "aspect_sin", "dist_to_road", "dist_to_village"),
  sample_name = "Woody"
)


overlap_balance_all <- bind_rows(forest_balance_overlap, woody_balance_overlap) %>%
  mutate(
    abs_smd = abs(smd),
    balance_flag = case_when(
      abs_smd < 0.1 ~ "Good",
      abs_smd < 0.2 ~ "Acceptable",
      TRUE ~ "Imbalanced"
    )
  )

print(overlap_balance_all)

write.csv(
  overlap_counts_all,
  "./plot/dai_hani_overlap_counts.csv",
  row.names = FALSE
)

write.csv(
  overlap_balance_all,
  "./plot/dai_hani_overlap_balance.csv",
  row.names = FALSE
)

overlap_maintext_table <- overlap_counts_all %>%
  tidyr::pivot_wider(
    names_from = c(stage, ethnicity),
    values_from = n
  ) %>%
  left_join(
    overlap_balance_all %>%
      group_by(sample, stage) %>%
      summarise(
        mean_abs_smd = mean(abs_smd, na.rm = TRUE),
        max_abs_smd = max(abs_smd, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = stage,
        values_from = c(mean_abs_smd, max_abs_smd)
      ),
    by = "sample"
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  ) %>%
  arrange(sample)

write.csv(
  overlap_maintext_table,
  "./plot/dai_hani_overlap_maintext.csv",
  row.names = FALSE
)

message("Saved overlap diagnostics:")
message("./plot/dai_hani_overlap_counts.csv")
message("./plot/dai_hani_overlap_balance.csv")
message("./plot/dai_hani_overlap_maintext.csv")





# ----------------------------------------------------------------
# 3. Nested models: does the ethnic SNS gap shrink after controls?
# ----------------------------------------------------------------
# Main logic:
# - If is_sns:dai shrinks sharply after geography controls or overlap restriction,
#   geography explains much of the Dai-Hani gap.
# - If is_sns:dai stays large and stable after strict restriction,
#   institutional difference becomes more plausible.

fit_nested_lm <- function(df, response, use_agri = TRUE) {
  base <- paste0(response, " ~ is_sns * dai_eth")
  terrain <- paste0(base, " + mean_elev + mean_slope")
  access <- paste0(terrain, " + dist_to_road + dist_to_village")

  socio <- if (use_agri) {
    paste0(access, " + log_income + log_popu + log_agri")
  } else {
    paste0(access, " + log_income + log_popu")
  }

  df_lm <- df %>%
    st_drop_geometry() %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      log_income = log(net_income_adjusted),
      log_popu = log(popu_density),
      log_agri = log(agri_land_pp)
    )

  list(
    M0 = lm(as.formula(base), data = df_lm),
    M1 = lm(as.formula(terrain), data = df_lm),
    M2 = lm(as.formula(access), data = df_lm),
    M3 = lm(as.formula(socio), data = df_lm)
  )
}

forest_nested <- fit_nested_lm(
  forest_overlap$overlap_df,
  "loss_1_345_intensity_14",
  use_agri = TRUE
)


woody_nested <- fit_nested_lm(
  woody_overlap$overlap_df,
  "loss_12_34_intensity_14",
  use_agri = TRUE
)

extract_term_path <- function(model_list, term_name) {
  bind_rows(lapply(names(model_list), function(mn) {
    sm <- summary(model_list[[mn]])$coef
    if (!(term_name %in% rownames(sm))) return(NULL)
    data.frame(
      model_step = mn,
      term = term_name,
      estimate = sm[term_name, "Estimate"],
      p_value = sm[term_name, grep("Pr\\(", colnames(sm))[1]],
      row.names = NULL
    )
  }))
}

forest_gap_path <- extract_term_path(forest_nested, "is_sns:dai_eth")
woody_gap_path  <- extract_term_path(woody_nested, "is_sns:dai_eth")

vif(forest_nested$M3)
vif(woody_nested$M3)

forest_is_sns_path <- extract_term_path(forest_nested, "is_sns")
woody_is_sns_path  <- extract_term_path(woody_nested, "is_sns")


forest_key_path <- bind_rows(forest_is_sns_path, forest_gap_path) %>%
  arrange(term, model_step)

woody_key_path <- bind_rows(woody_is_sns_path, woody_gap_path) %>%
  arrange(term, model_step)

write.csv(forest_key_path, "./plot/forest_key_terms_path.csv", row.names = FALSE)
write.csv(woody_key_path, "./plot/woody_key_terms_path.csv", row.names = FALSE)

# --------------------------------------------------------------------
# 4. Spatial version on overlap sample: if needed, re-run with lagsarlm
# --------------------------------------------------------------------

make_listw_knn <- function(sf_df, k = 5) {
  centroids <- st_centroid(sf_df)
  coords <- as.matrix(st_coordinates(centroids)[, 1:2])
  knn_obj <- knearneigh(coords, k = k)
  nb_obj <- knn2nb(knn_obj)
  nb2listw(nb_obj, style = "W")
}

fit_overlap_slm <- function(df, response, use_agri = TRUE) {
  listw_obj <- make_listw_knn(df, k = 5)

  df_fit <- df %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      log_income = log(net_income_adjusted),
      log_popu = log(popu_density),
      log_agri = log(agri_land_pp)
    )

  rhs <- if (use_agri) {
    paste(
      "is_sns * dai_eth + mean_elev + mean_slope + dist_to_road + dist_to_village +",
      "log_income + log_popu + log_agri"
    )
  } else {
    paste(
      "is_sns * dai_eth + mean_elev + mean_slope + dist_to_road + dist_to_village +",
      "log_income + log_popu"
    )
  }

  spatialreg::lagsarlm(
    as.formula(paste(response, "~", rhs)),
    data = df_fit,
    listw = listw_obj
  )
}

# Example:
forest_slm_overlap <- fit_overlap_slm(
   forest_overlap$overlap_df,
   "loss_1_345_intensity_14",
   use_agri = TRUE
 )
summary(forest_slm_overlap)

# ------------------------------------------------------------------
# 5. Geography moderation: does environment itself change SNS effect?
# ------------------------------------------------------------------
# This tests whether ethnic differences disappear once we allow
# sacred-site protection to vary by terrain/accessibility.

fit_geography_moderation <- function(df, response, use_agri = TRUE) {
  df_mod <- df %>%
    st_drop_geometry() %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      log_income = log(net_income_adjusted),
      log_popu = log(popu_density),
      log_agri = log(agri_land_pp)
    )

  rhs <- if (use_agri) {
    paste(
      "is_sns * dai_eth +",
      "is_sns:mean_elev + is_sns:dist_to_road +",
      "mean_elev + mean_slope + dist_to_road + dist_to_village +",
      "log_income + log_popu + log_agri"
    )
  } else {
    paste(
      "is_sns * dai_eth +",
      "is_sns:mean_elev + is_sns:dist_to_road +",
      "mean_elev + mean_slope + dist_to_road + dist_to_village +",
      "log_income + log_popu"
    )
  }

  lm(as.formula(paste(response, "~", rhs)), data = df_mod)
}

# Example:
woody_mod <- fit_geography_moderation(
   woody_overlap$overlap_df,
   "loss_12_34_intensity_14",
   use_agri = TRUE
 )
summary(woody_mod)


forest_mod <- fit_geography_moderation(
  forest_overlap$overlap_df,
  "loss_1_345_intensity_14",
  use_agri = TRUE
)
summary(forest_mod)

vif(forest_mod)
# ----------------------------------------------------------------------
# 6. Interpretation rule of thumb
# ----------------------------------------------------------------------
# A. If Dai-Hani true sites are geographically very different:
#    -> geography is a serious confounder and must be addressed directly
#
# B. If is_sns:dai_eth shrinks strongly after overlap restriction:
#    -> much of the ethnic gap is likely geography-driven
#
# C. If is_sns:dai_eth stays large after overlap + controls:
#    -> institutional interpretation becomes stronger
#
# D. If is_sns:geography terms are strong and ethnic gap weakens:
#    -> "institution × environment fit" may explain the pattern better
#
# Suggested next reporting sequence:
# 1. Show Dai-Hani true-site geography distributions
# 2. Show overlap/common-support restriction
# 3. Show nested-model coefficient path for is_sns:dai_eth
# 4. Show geography-moderation robustness
