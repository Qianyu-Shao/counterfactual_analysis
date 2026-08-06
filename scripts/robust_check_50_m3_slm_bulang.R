

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Locate the script and project root so the file can be run either from
# RStudio, the terminal, or another working directory.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
}

if (!is.null(script_file) && nzchar(script_file)) {
  script_file <- gsub("~\\+~", " ", script_file, fixed = FALSE)
}

project_root <- if (!is.null(script_file) && nzchar(script_file)) {
  normalizePath(
    file.path(dirname(script_file), ".."),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  normalizePath(
    "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis",
    winslash = "/",
    mustWork = TRUE
  )
}

setwd(project_root)

# Try common UTF-8 locale settings so Chinese file paths and shapefile fields
# are read consistently across macOS/R sessions.
invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

# Load the spatial, data-wrangling, table, and spatial-regression packages used
# throughout the repeated-draw workflow.
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(knitr)
  library(kableExtra)
  library(spdep)
  library(spatialreg)
})

# Define the repeated-draw range and reproducible seed rule. Each draw uses
# draw_seed = base_seed + draw_id, matching the alt-coding companion script.
start_draw <- 1
end_draw <- 50
n_draws <- end_draw - start_draw + 1
base_seed <- 20260509

# Define input script locations, output folders, and temporary raster folders.
# The main outputs are kept separate from the alternative-coding outputs.
scripts_dir <- file.path(project_root, "scripts")
output_dir <- file.path(project_root, "plot", "robust_check_50_m3_slm_bulang")
temp_root <- file.path(project_root, "plot", "terra_tmp_50_m3_slm_bulang")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temp_root, recursive = TRUE, showWarnings = FALSE)

coef_checkpoint_file <- file.path(output_dir, "draw_level_coefficients_checkpoint.csv")
counts_checkpoint_file <- file.path(output_dir, "draw_level_sample_counts_checkpoint.csv")
error_checkpoint_file <- file.path(output_dir, "draw_errors_checkpoint.csv")
stop_on_error <- TRUE

# Create a clean environment for each draw so objects created by sourced scripts
# do not accidentally carry over from the previous pseudo-site iteration.
reset_analysis_env <- function() {
  new.env(parent = globalenv())
}

# Source an R script into the draw-specific environment.
safe_source <- function(file, envir) {
  sys.source(file, envir = envir)
}

# Keep only Dai and Hani observations and add explicit ethnicity indicators used
# in overlap restriction and model formulas.
prepare_overlap_sample <- function(data) {
  data %>%
    filter(dai == 1 | hani == 1) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

# Count true and pseudo sites by ethnicity before and after the overlap
# restriction. These counts are later summarized across the 50 draws.
count_by_type <- function(df, sample_name, stage, draw_id) {
  counts <- df %>%
    st_drop_geometry() %>%
    mutate(
      ethnicity = ifelse(dai == 1, "Dai", "Hani"),
      site_type = ifelse(is_sns == 1, "True", "Pseudo")
    ) %>%
    count(ethnicity, site_type, name = "n")

  get_n <- function(eth, typ) {
    val <- counts$n[counts$ethnicity == eth & counts$site_type == typ]
    if (length(val) == 0) 0L else as.integer(val[1])
  }

  tibble(
    draw = draw_id,
    sample = sample_name,
    stage = stage,
    n_dai_true = get_n("Dai", "True"),
    n_hani_true = get_n("Hani", "True"),
    n_dai_pseudo = get_n("Dai", "Pseudo"),
    n_hani_pseudo = get_n("Hani", "Pseudo")
  )
}

# Restrict Dai and Hani observations to their common geographic support. A
# propensity score for being Dai is estimated using terrain and accessibility
# covariates, then observations outside the overlapping Dai-Hani score range are
# removed.
build_ethnicity_overlap <- function(df) {
  base_df <- df %>%
    st_drop_geometry() %>%
    mutate(
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180)
    )

  ps_fit <- glm(
    dai_eth ~ mean_elev + mean_slope + aspect_cos + aspect_sin +
      dist_to_road + dist_to_village,
    family = binomial(),
    data = base_df
  )

  base_df$ps_eth <- predict(ps_fit, type = "response")

  dai_range <- range(base_df$ps_eth[base_df$dai_eth == 1], na.rm = TRUE)
  hani_range <- range(base_df$ps_eth[base_df$dai_eth == 0], na.rm = TRUE)

  lower <- max(dai_range[1], hani_range[1])
  upper <- min(dai_range[2], hani_range[2])

  overlap_ids <- base_df %>%
    filter(ps_eth >= lower, ps_eth <= upper) %>%
    pull(ID)

  df %>%
    filter(ID %in% overlap_ids) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

# Build a row-standardized k-nearest-neighbor spatial weights object from site
# centroids. This is the spatial W matrix used by the spatial lag model.
make_listw_knn <- function(sf_df, k = 5) {
  centroids <- st_centroid(sf_df)
  coords <- as.matrix(st_coordinates(centroids)[, 1:2])
  nb_obj <- knn2nb(knearneigh(coords, k = k))
  nb2listw(nb_obj, style = "W")
}

# Fit one M3 spatial lag model for a given outcome. The formula includes the
# sacred-site indicator, ethnicity, their interaction, terrain/accessibility
# controls, circular aspect controls, and village-level socioeconomic controls.
fit_m3_slm <- function(df, response, ethnicity_var, sample_name, outcome_name, draw_id, k = 5) {
  df_fit <- df %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180),
      log_income = log(net_income_adjusted),
      log_popu = log(popu_density),
      log_agri = log(agri_land_pp)
    )

  rhs <- paste(
    "is_sns *", ethnicity_var, "+",
    "mean_elev + mean_slope + aspect_cos + aspect_sin +",
    "dist_to_road + dist_to_village +",
    "log_income + log_popu + log_agri"
  )

  formula_obj <- as.formula(paste(response, "~", rhs))
  listw_obj <- make_listw_knn(df, k = k)

  slm_fit <- spatialreg::lagsarlm(formula_obj, data = df_fit, listw = listw_obj)

  list(
    sample = sample_name,
    outcome = outcome_name,
    draw = draw_id,
    slm = slm_fit
  )
}

# Convert the spatial lag model summary into a tidy coefficient table and add
# rho, the spatial autoregressive parameter, as an additional term.
extract_all_coefs <- function(model, sample_name, outcome_name, draw_id) {
  coef_mat <- as.data.frame(summary(model)$Coef)
  coef_mat$term <- rownames(coef_mat)
  rownames(coef_mat) <- NULL

  est_col <- grep("Estimate", names(coef_mat), value = TRUE)[1]
  p_col <- grep("Pr\\(", names(coef_mat), value = TRUE)[1]

  out <- coef_mat %>%
    transmute(
      draw = draw_id,
      sample = sample_name,
      outcome = outcome_name,
      term = term,
      estimate = as.numeric(.data[[est_col]]),
      p_value = as.numeric(.data[[p_col]])
    )

  rho_term <- if (!is.null(model$rho)) {
    tibble(
      draw = draw_id,
      sample = sample_name,
      outcome = outcome_name,
      term = "rho",
      estimate = as.numeric(model$rho),
      p_value = NA_real_
    )
  } else {
    NULL
  }

  bind_rows(out, rho_term)
}

# Save intermediate coefficient/count/error files after every draw. This makes
# the long 50-draw run resumable/inspectable if a later draw fails.
write_checkpoint <- function(coef_results, count_results, draw_errors) {
  coef_df <- bind_rows(coef_results)
  count_df <- bind_rows(count_results)
  error_df <- bind_rows(draw_errors)

  if (nrow(coef_df) > 0) write_csv(coef_df, coef_checkpoint_file)
  if (nrow(count_df) > 0) write_csv(count_df, counts_checkpoint_file)
  if (nrow(error_df) > 0) write_csv(error_df, error_checkpoint_file)
}

# Run the full pipeline for one pseudo-site draw: source data-preparation
# scripts, create matched samples, apply Dai-Hani overlap restriction, count
# samples, fit the four outcome-specific SLMs, and return tidy outputs.
run_one_draw <- function(draw_id) {
  draw_env <- reset_analysis_env()
  draw_env$.Random.seed <- NULL
  draw_seed <- base_seed + draw_id
  set.seed(draw_seed)

  if (requireNamespace("terra", quietly = TRUE)) {
    terra::terraOptions(tempdir = temp_root)
  }
  if (requireNamespace("raster", quietly = TRUE)) {
    raster::rasterOptions(tmpdir = temp_root)
  }

  safe_source(file.path(scripts_dir, "function.R"), draw_env)
  safe_source(file.path(scripts_dir, "bulang_clean.R"), draw_env)
  safe_source(file.path(scripts_dir, "within_ethnicity_match.R"), draw_env)
  safe_source(file.path(scripts_dir, "extract_encroachment_after_match.R"), draw_env)

  forest_before <- prepare_overlap_sample(draw_env$matched_forest_dai_hani)
  woody_before <- prepare_overlap_sample(draw_env$matched_woody_dai_hani)

  forest_after <- build_ethnicity_overlap(forest_before)
  woody_after <- build_ethnicity_overlap(woody_before)

  count_df <- bind_rows(
    count_by_type(forest_before, "Forest", "Before overlap", draw_id),
    count_by_type(forest_after, "Forest", "After overlap", draw_id),
    count_by_type(woody_before, "Woody", "Before overlap", draw_id),
    count_by_type(woody_after, "Woody", "After overlap", draw_id)
  ) %>%
    mutate(draw_seed = draw_seed)

  forest_cashcrop <- fit_m3_slm(
    forest_after,
    response = "loss_1_34_intensity_14",
    ethnicity_var = "hani_eth",
    sample_name = "Forest",
    outcome_name = "forest_to_cashcrop",
    draw_id = draw_id,
    k = 5
  )

  forest_farmland <- fit_m3_slm(
    forest_after,
    response = "loss_1_345_intensity_14",
    ethnicity_var = "hani_eth",
    sample_name = "Forest",
    outcome_name = "forest_to_farmland",
    draw_id = draw_id,
    k = 5
  )

  woody_cashcrop <- fit_m3_slm(
    woody_after,
    response = "loss_12_34_intensity_14",
    ethnicity_var = "dai_eth",
    sample_name = "Woody",
    outcome_name = "woody_to_cashcrop",
    draw_id = draw_id,
    k = 5
  )

  woody_farmland <- fit_m3_slm(
    woody_after,
    response = "loss_12_345_intensity_14",
    ethnicity_var = "dai_eth",
    sample_name = "Woody",
    outcome_name = "woody_to_farmland",
    draw_id = draw_id,
    k = 5
  )

  coef_df <- bind_rows(
    extract_all_coefs(forest_cashcrop$slm, "Forest", "forest_to_cashcrop", draw_id),
    extract_all_coefs(forest_farmland$slm, "Forest", "forest_to_farmland", draw_id),
    extract_all_coefs(woody_cashcrop$slm, "Woody", "woody_to_cashcrop", draw_id),
    extract_all_coefs(woody_farmland$slm, "Woody", "woody_to_farmland", draw_id)
  ) %>%
    mutate(draw_seed = draw_seed)

  list(coefs = coef_df, counts = count_df)
}

# Summarize coefficient stability across repeated draws, including sign
# consistency and the share of draws with p < 0.05 or p < 0.10.
summarize_coef_stability <- function(df) {
  df %>%
    group_by(sample, outcome, term) %>%
    summarise(
      n_success = sum(!is.na(estimate)),
      mean_estimate = mean(estimate, na.rm = TRUE),
      positive_share = mean(estimate > 0, na.rm = TRUE),
      negative_share = mean(estimate < 0, na.rm = TRUE),
      p_lt_005_rate = mean(p_value < 0.05, na.rm = TRUE),
      p_lt_010_rate = mean(p_value < 0.10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        c(mean_estimate, positive_share, negative_share, p_lt_005_rate, p_lt_010_rate),
        ~ round(.x, 3)
      )
    ) %>%
    arrange(sample, outcome, term)
}

# Summarize how sample sizes vary before and after overlap restriction across
# the 50 repeated pseudo-site draws.
summarize_counts <- function(df) {
  df %>%
    group_by(sample, stage) %>%
    summarise(
      mean_n_dai_true = mean(n_dai_true, na.rm = TRUE),
      mean_n_hani_true = mean(n_hani_true, na.rm = TRUE),
      mean_n_dai_pseudo = mean(n_dai_pseudo, na.rm = TRUE),
      mean_n_hani_pseudo = mean(n_hani_pseudo, na.rm = TRUE),
      min_n_dai_true = min(n_dai_true, na.rm = TRUE),
      min_n_hani_true = min(n_hani_true, na.rm = TRUE),
      min_n_dai_pseudo = min(n_dai_pseudo, na.rm = TRUE),
      min_n_hani_pseudo = min(n_hani_pseudo, na.rm = TRUE),
      max_n_dai_true = max(n_dai_true, na.rm = TRUE),
      max_n_hani_true = max(n_hani_true, na.rm = TRUE),
      max_n_dai_pseudo = max(n_dai_pseudo, na.rm = TRUE),
      max_n_hani_pseudo = max(n_hani_pseudo, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(-c(sample, stage), ~ round(.x, 1)))
}

# Format summary tables for quick HTML inspection.
format_html_table <- function(df, caption) {
  knitr::kable(
    df,
    format = "html",
    caption = caption
  ) %>%
    kableExtra::kable_styling(full_width = FALSE)
}

message(sprintf("Starting 50-draw M3 SLM robustness check with %s draws...", n_draws))

# Allocate list containers for draw-level outputs before entering the loop.
coef_results <- vector("list", n_draws)
count_results <- vector("list", n_draws)
draw_errors <- vector("list", n_draws)

# Clear temporary raster files from previous runs before starting a new 50-draw
# sequence.
unlink(file.path(temp_root, "*"), recursive = TRUE, force = TRUE)

# Main repeated-draw loop. Each iteration uses its own seed, writes checkpoints,
# and either stops on failure or stores successful draw-level outputs.
for (idx in seq_len(n_draws)) {
  draw_id <- start_draw + idx - 1
  message(sprintf("Running draw %s / %s (global draw %s)", idx, n_draws, draw_id))

  result <- tryCatch(
    run_one_draw(draw_id),
    error = function(e) {
      draw_errors[[idx]] <<- tibble(draw = draw_id, error = conditionMessage(e))
      NULL
    }
  )

  if (is.null(result)) {
    write_checkpoint(coef_results, count_results, draw_errors)
    if (isTRUE(stop_on_error)) {
      stop(sprintf("Draw %s failed: %s", draw_id, draw_errors[[idx]]$error))
    }
  } else {
    coef_results[[idx]] <- result$coefs
    count_results[[idx]] <- result$counts
    write_checkpoint(coef_results, count_results, draw_errors)
  }
}

# Combine all draw-level results and build across-draw summary tables.
coef_df <- bind_rows(coef_results)
count_df <- bind_rows(count_results)
error_df <- bind_rows(draw_errors)

coef_summary_df <- summarize_coef_stability(coef_df)
count_summary_df <- summarize_counts(count_df)

# Write machine-readable CSV outputs for downstream figures and manuscript
# tables.
write_csv(coef_df, file.path(output_dir, "draw_level_coefficients.csv"))
write_csv(count_df, file.path(output_dir, "draw_level_sample_counts.csv"))
write_csv(coef_summary_df, file.path(output_dir, "coef_summary.csv"))
write_csv(count_summary_df, file.path(output_dir, "sample_count_summary.csv"))

if (nrow(error_df) > 0) {
  write_csv(error_df, file.path(output_dir, "draw_errors.csv"))
}

# Write HTML versions of the two main summaries for quick visual checking.
coef_html <- format_html_table(
  coef_summary_df,
  "M3 SLM coefficient stability across 50 repeated pseudo-site draws"
)

count_html <- format_html_table(
  count_summary_df,
  "Sample counts before and after overlap across 50 repeated pseudo-site draws"
)

writeLines(as.character(coef_html), file.path(output_dir, "coef_summary.html"))
writeLines(as.character(count_html), file.path(output_dir, "sample_count_summary.html"))

# Save the complete analysis object, including draw-level outputs, summaries,
# errors, and configuration metadata.
saveRDS(
  list(
    coefficients = coef_df,
    sample_counts = count_df,
    coef_summary = coef_summary_df,
    count_summary = count_summary_df,
    draw_errors = error_df,
    config = list(
      start_draw = start_draw,
      end_draw = end_draw,
      n_draws = n_draws,
      base_seed = base_seed,
      seed_rule = "draw_seed = base_seed + draw_id",
      source_script = "bulang_clean.R",
      analysis_target = "M3 overlap sample SLM only, four outcomes"
    )
  ),
  file.path(output_dir, "robust_m3_slm_50.rds")
)

message("Saved outputs to:")
message(output_dir)
