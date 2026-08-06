`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

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
  normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    "/Users/shaoyuchen/Desktop/thesis/pseudo sites/counterfactual_analysis",
    winslash = "/",
    mustWork = TRUE
  )
}

setwd(project_root)

invisible(try(Sys.setlocale("LC_CTYPE", "zh_CN.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE))
invisible(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE))

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(tibble)
})

start_draw <- 1
end_draw <- 50
base_seed <- 20260509

scripts_dir <- file.path(project_root, "scripts")
output_root <- file.path(project_root, "plot", "matched_overlap_inputs_50_draws")
forest_dir <- file.path(output_root, "forest")
woody_dir <- file.path(output_root, "woody")
temp_root <- file.path(project_root, "plot", "terra_tmp_export_50_matched_overlap")

dir.create(forest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(woody_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temp_root, recursive = TRUE, showWarnings = FALSE)

reset_analysis_env <- function() {
  new.env(parent = globalenv())
}

safe_source <- function(file, envir) {
  sys.source(file, envir = envir)
}

prepare_overlap_sample <- function(data) {
  data %>%
    filter(dai == 1 | hani == 1) %>%
    mutate(
      dai_eth = ifelse(dai == 1, 1, 0),
      hani_eth = ifelse(hani == 1, 1, 0),
      ethnicity = ifelse(dai == 1, "Dai", "Hani")
    )
}

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

add_regression_variables <- function(df, draw_id, draw_seed, sample_name) {
  df %>%
    mutate(
      draw = draw_id,
      draw_seed = draw_seed,
      sample = sample_name,
      aspect_cos = cos(mode_aspect * pi / 180),
      aspect_sin = sin(mode_aspect * pi / 180),
      log_income = log(net_income_adjusted),
      log_popu = log(popu_density),
      log_agri = log(agri_land_pp),
      geometry_wkt = sf::st_as_text(sf::st_geometry(.))
    ) %>%
    st_drop_geometry() %>%
    relocate(draw, draw_seed, sample, .before = 1)
}

export_one_draw <- function(draw_id) {
  draw_seed <- base_seed + draw_id
  set.seed(draw_seed)

  if (requireNamespace("terra", quietly = TRUE)) {
    terra::terraOptions(tempdir = temp_root)
  }
  if (requireNamespace("raster", quietly = TRUE)) {
    raster::rasterOptions(tmpdir = temp_root)
  }

  draw_env <- reset_analysis_env()
  safe_source(file.path(scripts_dir, "function.R"), draw_env)
  safe_source(file.path(scripts_dir, "bulang_clean.R"), draw_env)
  safe_source(file.path(scripts_dir, "within_ethnicity_match.R"), draw_env)
  safe_source(file.path(scripts_dir, "extract_encroachment_after_match.R"), draw_env)

  forest_after <- draw_env$matched_forest_dai_hani %>%
    prepare_overlap_sample() %>%
    build_ethnicity_overlap() %>%
    add_regression_variables(draw_id, draw_seed, "Forest")

  woody_after <- draw_env$matched_woody_dai_hani %>%
    prepare_overlap_sample() %>%
    build_ethnicity_overlap() %>%
    add_regression_variables(draw_id, draw_seed, "Woody")

  write_csv(
    forest_after,
    file.path(forest_dir, sprintf("draw_%02d_forest_matched_overlap.csv", draw_id))
  )

  write_csv(
    woody_after,
    file.path(woody_dir, sprintf("draw_%02d_woody_matched_overlap.csv", draw_id))
  )

  tibble(
    draw = draw_id,
    draw_seed = draw_seed,
    forest_n = nrow(forest_after),
    woody_n = nrow(woody_after),
    forest_file = file.path(forest_dir, sprintf("draw_%02d_forest_matched_overlap.csv", draw_id)),
    woody_file = file.path(woody_dir, sprintf("draw_%02d_woody_matched_overlap.csv", draw_id))
  )
}

message("Exporting matched + overlap regression input datasets...")
unlink(file.path(temp_root, "*"), recursive = TRUE, force = TRUE)

manifest <- lapply(start_draw:end_draw, function(draw_id) {
  message(sprintf("Exporting draw %02d / %02d", draw_id, end_draw))
  export_one_draw(draw_id)
}) %>%
  bind_rows()

write_csv(manifest, file.path(output_root, "matched_overlap_input_manifest.csv"))

message("Saved 50 forest CSVs and 50 woody CSVs to: ", output_root)
