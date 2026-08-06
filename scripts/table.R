extract_model_info <- function(model, name) {
  model_type <- if (!is.null(model$rho)) {
    "SLM"
  } else if (!is.null(model$lambda)) {
    "SEM"
  } else {
    "Unknown"
  }
  
  spatial_value <- if (!is.null(model$rho)) {
    model$rho
  } else if (!is.null(model$lambda)) {
    model$lambda
  } else {
    NA_real_
  }
  
  list(
    model = name,
    coef = summary(model)$Coef,
    n = length(residuals(model)),
    aic = AIC(model),
    logLik = as.numeric(logLik(model)),
    model_type = model_type,
    spatial_name = "Spatial",
    spatial_value = spatial_value
  )
}

model_info <- list(
  extract_model_info(slm_fit_1_3, "forest to cashcrop"),
  extract_model_info(slm_fit_1_5, "forest to farmland"),
  extract_model_info(slm_fit_12_3, "woody land to cashcrop"),
  extract_model_info(slm_fit_12_5, "woody land to farmland")
)

library(dplyr)
library(tidyr)
library(tibble)
library(knitr)
library(kableExtra)
library(gridExtra)
library(grid)

term_labels <- c(
  "(Intercept)" = "Intercept",
  "is_sns" = "Is SNS",
  "dai" = "Dai",
  "is_sns:dai" = "SNS x Dai",
  "log_income" = "Log Income",
  "log_popu" = "Log Population",
  "log_agri" = "Log Agricultural Land"
)

add_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ ""
  )
}

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0 || is.null(x) || is.na(x[1])) {
    return("\u2014")
  }
  sprintf(paste0("%.", digits, "f"), as.numeric(x[1]))
}

extract_coef_rows <- function(info, term_labels = NULL) {
  coef_mat <- as.data.frame(info$coef)
  coef_mat$term_raw <- rownames(coef_mat)
  rownames(coef_mat) <- NULL
  
  p_col <- grep("Pr\\(", names(coef_mat), value = TRUE)[1]
  est_col <- grep("Estimate", names(coef_mat), value = TRUE)[1]
  
  out <- coef_mat %>%
    mutate(
      term = if (!is.null(term_labels)) {
        ifelse(term_raw %in% names(term_labels), term_labels[term_raw], term_raw)
      } else {
        term_raw
      },
      order_id = row_number()
    ) %>%
    transmute(
      order_id,
      model = info$model,
      coef_term = term,
      p_term = paste0(term, "_p"),
      coef_value = sprintf("%.3f%s", .data[[est_col]], add_stars(.data[[p_col]])),
      p_value = sprintf("(%.3f)", .data[[p_col]])
    )
  
  bind_rows(
    out %>% transmute(order_id = order_id * 2 - 1, term = coef_term, model, value = coef_value),
    out %>% transmute(order_id = order_id * 2, term = p_term, model, value = p_value)
  )
}

extract_stat_rows <- function(info, start_id) {
  tibble(
    order_id = start_id + 1:5,
    term = c("N", "AIC", "LogLik", "Model Type", info$spatial_name),
    model = info$model,
    value = c(
      as.character(info$n),
      fmt_num(info$aic),
      fmt_num(info$logLik),
      info$model_type,
      fmt_num(info$spatial_value)
    )
  )
}

coef_long <- bind_rows(lapply(model_info, extract_coef_rows, term_labels = term_labels))
max_id <- max(coef_long$order_id)
stats_long <- bind_rows(lapply(model_info, extract_stat_rows, start_id = max_id))

result_wide <- bind_rows(coef_long, stats_long) %>%
  arrange(order_id) %>%
  select(order_id, term, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  select(-order_id) %>%
  mutate(
    term = gsub("_p$", "", term),
    term = ifelse(duplicated(term), "", term)
  )

table_html <- kable(
  result_wide,
  format = "html",
  escape = FALSE,
  col.names = c("Term", vapply(model_info, `[[`, "", "model")),
  caption = "Coefficients and model statistics"
) %>%
  kable_styling(full_width = FALSE) %>%
  footnote(
    general = "Entries report coefficients with significance stars; p-values are shown in parentheses on the row below. Spatial reports rho for SLM and lambda for SEM. Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1.",
    general_title = "Note: "
  )

writeLines(as.character(table_html), "./plot/slm_summary_table.html")

png(
  filename = "./plot/slm_summary_table.png",
  width = 2600,
  height = 1900,
  res = 220
)
grid.newpage()
grid.draw(
  tableGrob(
    result_wide,
    rows = NULL,
    theme = ttheme_minimal(
      base_size = 12,
      core = list(fg_params = list(hjust = c(0, rep(0.5, ncol(result_wide) - 1)))),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
)
grid.text(
  "Note: Coefficients with significance stars; p-values are shown in parentheses on the row below. Spatial reports rho for SLM and lambda for SEM. *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1.",
  x = unit(0.01, "npc"),
  y = unit(0.02, "npc"),
  just = c("left", "bottom"),
  gp = gpar(fontsize = 9)
)
dev.off()

message("Saved:")
message("./plot/slm_summary_table.html")
