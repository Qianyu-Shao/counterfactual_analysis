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

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
})

source_dirs <- tibble(
  source_version = c("bulang", "bulang_altcoding"),
  description = c(
    "Main Dai-Hani coding used for ethnic-specific is_sns effects",
    "Alternative Dai-Hani coding used for is_sns:dai interaction effects"
  ),
  source_dir = file.path(
    project_root,
    "plot",
    c("robust_check_50_m3_slm_bulang", "robust_check_50_m3_slm_bulang_altcoding")
  )
)

collection_dir <- file.path(project_root, "plot", "repeated_spatiallag_50_draw_data")
raw_dir <- file.path(collection_dir, "source_outputs")
dir.create(collection_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "draw_level_coefficients.csv",
  "draw_level_sample_counts.csv",
  "sample_count_summary.csv",
  "coef_summary.csv",
  "robust_m3_slm_50.rds"
)

optional_files <- c(
  "draw_level_coefficients_checkpoint.csv",
  "draw_level_sample_counts_checkpoint.csv",
  "draw_errors_checkpoint.csv"
)

copy_source_outputs <- function(source_version, source_dir) {
  target_dir <- file.path(raw_dir, source_version)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  files_to_copy <- c(required_files, optional_files)
  source_paths <- file.path(source_dir, files_to_copy)
  keep <- file.exists(source_paths)

  file.copy(
    from = source_paths[keep],
    to = file.path(target_dir, files_to_copy[keep]),
    overwrite = TRUE
  )

  tibble(
    source_version = source_version,
    file = files_to_copy,
    source_path = source_paths,
    copied_to = file.path(target_dir, files_to_copy),
    available = keep
  )
}

copy_manifest <- pmap_dfr(
  source_dirs %>% select(source_version, source_dir),
  copy_source_outputs
)

read_tagged_csv <- function(source_version, source_dir, file_name) {
  path <- file.path(source_dir, file_name)
  if (!file.exists(path)) {
    return(tibble())
  }

  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      source_version = source_version,
      .before = 1
    )
}

combined_coefficients <- pmap_dfr(
  source_dirs %>% select(source_version, source_dir),
  ~ read_tagged_csv(..1, ..2, "draw_level_coefficients.csv")
)

combined_sample_counts <- pmap_dfr(
  source_dirs %>% select(source_version, source_dir),
  ~ read_tagged_csv(..1, ..2, "draw_level_sample_counts.csv")
)

combined_coef_summary <- pmap_dfr(
  source_dirs %>% select(source_version, source_dir),
  ~ read_tagged_csv(..1, ..2, "coef_summary.csv")
)

combined_sample_count_summary <- pmap_dfr(
  source_dirs %>% select(source_version, source_dir),
  ~ read_tagged_csv(..1, ..2, "sample_count_summary.csv")
)

draw_metadata <- combined_coefficients %>%
  distinct(source_version, draw, draw_seed) %>%
  arrange(source_version, draw)

write_csv(source_dirs, file.path(collection_dir, "source_versions.csv"))
write_csv(copy_manifest, file.path(collection_dir, "collection_manifest.csv"))
write_csv(draw_metadata, file.path(collection_dir, "draw_seed_metadata.csv"))
write_csv(combined_coefficients, file.path(collection_dir, "combined_draw_level_coefficients.csv"))
write_csv(combined_sample_counts, file.path(collection_dir, "combined_draw_level_sample_counts.csv"))
write_csv(combined_coef_summary, file.path(collection_dir, "combined_coef_summary.csv"))
write_csv(combined_sample_count_summary, file.path(collection_dir, "combined_sample_count_summary.csv"))

message("Collected repeated spatial lag outputs in: ", collection_dir)
