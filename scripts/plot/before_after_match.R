library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

forest_before <- final_buffer_500_12345_forest %>%
  sf::st_drop_geometry() %>%
  mutate(
    ethnicity = case_when(
      dai == 1 ~ "Dai",
      hani == 1 ~ "Hani",
      TRUE ~ "Bulang"
    )
  ) %>%
  filter(ethnicity %in% c("Dai", "Hani", "Bulang")) %>%
  mutate(sample = "Before matching")

forest_after <- bind_rows(
  matched_forest_dai %>% sf::st_drop_geometry(),
  matched_forest_hani %>% sf::st_drop_geometry(),
  matched_forest_bulang %>% sf::st_drop_geometry()
) %>%
  mutate(sample = "After matching")

forest_vars <- c(
  "mean_elev",
  "mean_slope",
  "dist_to_road",
  "mode_aspect",
  "forest_area_02",
  "dist_to_village"
)

forest_plot_data <- bind_rows(forest_before, forest_after) %>%
  select(all_of(forest_vars), is_sns, sample) %>%
  pivot_longer(
    cols = all_of(forest_vars),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value), !is.na(is_sns)) %>%
  mutate(
    is_sns = factor(is_sns, levels = c(0, 1), labels = c("Pseudo-SNS", "True SNS")),
    variable = case_when(
      variable == "mean_elev" ~ "Elevation",
      variable == "mean_slope" ~ "Slope",
      variable == "dist_to_road" ~ "Distance to Road",
      variable == "mode_aspect" ~ "Aspect",
      variable == "forest_area_02" ~ "Initial Forest Area",
      variable == "dist_to_village" ~ "Distance to Village",
      TRUE ~ variable
    ),
    value_plot = case_when(
      variable %in% c("Distance to Road", "Distance to Village", "Initial Forest Area") ~ log1p(value),
      TRUE ~ value
    )
  )

p_forest <- ggplot(forest_plot_data, aes(x = value_plot, color = is_sns, fill = is_sns)) +
  geom_density(alpha = 0.2, adjust = 1.1) +
  facet_wrap(~ sample + variable, scales = "free", ncol = 3) +
  scale_color_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
  scale_fill_manual(values = c("Pseudo-SNS" = "#D95F02", "True SNS" = "#1B9E77")) +
  labs(
    title = "Within-Ethnicity Matching: Forest Sample",
    subtitle = "Combined before/after balance across Dai, Hani, and Bulang",
    x = NULL,
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  "./plot/within_ethnicity_before_after_forest.png",
  p_forest,
  width = 14,
  height = 10,
  dpi = 300
)
