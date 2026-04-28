# --------------------------------------------------------------------------------
# MGTWR: Comparative Analysis of Spatio-Temporal Calibration Algorithms
# --------------------------------------------------------------------------------

rm(list=ls())
gc()

# --------------------------------------------------------------------------------
# PACKAGE INSTALLATION & SETUP
# --------------------------------------------------------------------------------
if (!require(pacman)) install.packages("pacman")
if (!require(pacman)) install.packages("pacman")
pacman::p_load(                             
  dplyr, tidyr, purrr, stringr, readr,       # Data manipulation & text
  sf, spdep, geobr, censobr, classInt,       # Spatial data & weights
  ggplot2, ggspatial, ggridges, patchwork,   # Visualization & mapping
  ggdist, ggrepel, corrplot,                 # Plot geoms & stats
  gtsummary, gt, knitr, kableExtra,          # Table generation
  here, RANN, MASS                           # Utilities & math operations
)

# Output directories
dir.create(here::here("outputs", "simulated"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs", "monte_carlo"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs", "real"), recursive = TRUE, showWarnings = FALSE)

# Load functions
source(here::here("script", "Functions_code.R")) 

data_path <- here::here("data", "mgtwr_panel.rds")

map_theme_addons <- ggplot2::theme(
  axis.text = ggplot2::element_blank(), 
  axis.ticks = ggplot2::element_blank(), 
  panel.grid = ggplot2::element_blank()
)

# --------------------------------------------------------------------------------
# THE MONTE CARLO SIMULATION
# --------------------------------------------------------------------------------
N_sim <- 300
grid_sizes <- c(10, 15, 20) 

results_bias <- list()
results_metrics <- list()

formula_sim <- Y ~ X1 + X2

set.seed(42)
output_path <- here::here("outputs", "monte_carlo", "MC_Comparative_Results.rds")

if (file.exists(output_path)) {
  saved_results <- readRDS(output_path)
  all_bias <- saved_results$bias
  all_metrics <- saved_results$metrics
  
} else {
  message("No existing results found. Starting simulations...")
  
  results_metrics <- list()
  results_bias <- list()
  
  for (g_size in grid_sizes) {
    message(sprintf("\n------ Testing Grid Size: %dx%d (N = %d) ------", g_size, g_size, g_size^2 * 3))
    
    for (sim in 1:N_sim) {
      if (sim %% 100 == 0) message(sprintf(" >Running Simulation %d / %d", sim, N_sim))
      
      sim_data <- generate_dgp(seed = 1000 + sim, grid_size = g_size)
      coords_eval <- sf::st_coordinates(sf::st_centroid(sim_data))
      n_obs <- nrow(sim_data)
      
      # Models
      mod_2scall <- run_mgtwr_2scall(formula_sim, sim_data, time_var = "time") 
      mod_tds    <- run_mgtwr_tds(formula_sim, sim_data, time_var = "time")    
      mod_f2017  <- run_mgtwr_fotheringham(formula_sim, sim_data, time_var = "time") 
      
      # PIT for calibration metrics
      pit_2scall <- compute_continuous_pit(mod_2scall$residual)
      pit_tds    <- compute_continuous_pit(mod_tds$residual)
      pit_f2017  <- compute_continuous_pit(mod_f2017$residual)
      
      # Metrics
      metrics_2scall <- calc_model_metrics(sim_data$Y, mod_2scall$yhat, coords_eval, k = 5, exec_time = mod_2scall$exec_time)
      metrics_tds    <- calc_model_metrics(sim_data$Y, mod_tds$yhat, coords_eval, k = 5, exec_time = mod_tds$exec_time)
      metrics_f2017  <- calc_model_metrics(sim_data$Y, mod_f2017$yhat, coords_eval, k = 5, exec_time = mod_f2017$exec_time)
      
      # Metrics
      metrics_iter <- tibble::tibble(
        Sim_ID    = sim,
        Grid_Size = g_size,
        N_Obs     = n_obs,
        Model     = c("2SCALL", "TDS", "F2017"),
        RMSE      = c(metrics_2scall$RMSE, metrics_tds$RMSE, metrics_f2017$RMSE),
        Time_sec  = c(metrics_2scall$Time_sec, metrics_tds$Time_sec, metrics_f2017$Time_sec),
        U2        = c(watson_U2(pit_2scall)$U2, watson_U2(pit_tds)$U2, watson_U2(pit_f2017)$U2),
        S_star    = c(neyman_ledwina(pit_2scall)$Sstar, neyman_ledwina(pit_tds)$Sstar, neyman_ledwina(pit_f2017)$Sstar),
        TB2       = c(tb2_entropy_proxy(pit_2scall)$TB2, tb2_entropy_proxy(pit_tds)$TB2, tb2_entropy_proxy(pit_f2017)$TB2)
      )
      results_metrics[[length(results_metrics) + 1]] <- metrics_iter
      
      # Bias
      df_metrics <- sim_data |> sf::st_drop_geometry() |> 
        dplyr::select(tract_id, time, True_Beta_Intercept, True_Beta_X1, True_Beta_X2) |> 
        dplyr::mutate(
          Sim_ID    = sim,
          Grid_Size = g_size,
          N_Obs     = n_obs,
          
          Bias_Int_2SCALL = mod_2scall$coefs[,1] - True_Beta_Intercept,
          Bias_X1_2SCALL  = mod_2scall$coefs[,2] - True_Beta_X1,
          Bias_X2_2SCALL  = mod_2scall$coefs[,3] - True_Beta_X2,
          
          Bias_Int_TDS    = mod_tds$coefs[,1] - True_Beta_Intercept,
          Bias_X1_TDS     = mod_tds$coefs[,2] - True_Beta_X1,
          Bias_X2_TDS     = mod_tds$coefs[,3] - True_Beta_X2,
          
          Bias_Int_F2017  = mod_f2017$coefs[,1] - True_Beta_Intercept,
          Bias_X1_F2017   = mod_f2017$coefs[,2] - True_Beta_X1,
          Bias_X2_F2017   = mod_f2017$coefs[,3] - True_Beta_X2
        )
      results_bias[[length(results_bias) + 1]] <- df_metrics
    }
  }
  
  all_bias <- dplyr::bind_rows(results_bias)
  all_metrics <- dplyr::bind_rows(results_metrics)
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(list(bias = all_bias, metrics = all_metrics), output_path)
}

# --------------------------------------------------------------------------------
# EVALUATION PLOTTING (MONTE CARLO)
# --------------------------------------------------------------------------------
df_plot_bias <- all_bias |> 
  dplyr::select(N_Obs, starts_with("Bias_")) |> 
  tidyr::pivot_longer(-N_Obs, names_to = c("Parameter", "Model"), names_pattern = "Bias_(.*)_(.*)", values_to = "Bias") |>
  dplyr::mutate(
    Parameter = factor(Parameter, levels = c("Int", "X1", "X2"), labels = c("Intercept", "X1", "X2")),
    N_Label   = factor(paste0("N = ", N_Obs), levels = paste0("N = ", sort(unique(N_Obs)))),
    Model=recode(
      Model, "F2017"="BF"
      )
  )

# PARAMETER BIAS
#-------------------------------------------------------------------------------
p_box_scale <- ggplot(df_plot_bias, aes(x = Model, y = Bias, color = Model, fill = Model)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 1) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black", linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
  geom_hline(yintercept = 0, color = "#e74c3c", linetype = "dashed", linewidth = 1) +
  facet_grid(Parameter ~ N_Label, scales = "free_y") +
  scale_fill_viridis_d(option = "mako", begin = 0.1, end = 0.8) +
  scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.8) +
  theme_n() +
  labs(title = "",
       y = "Estimation Bias (Estimated - True)", x = "Calibration Algorithm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

ggsave(here::here("outputs", "monte_carlo", "Parameter Bias Recovery Across Sample Sizes.pdf"), plot = p_box_scale, width = 14, height = 8, dpi = 300, bg = "white")

#CONFIDENCE INTERVALS
#-------------------------------------------------------------------------------
df_ci <- df_plot_bias |>
  dplyr::group_by(Model, Parameter, N_Label) |>
  dplyr::summarise(
    Mean_Bias = mean(Bias, na.rm = TRUE),
    SD_Bias = sd(Bias, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    SE = SD_Bias / sqrt(N),
    CI_Lower = Mean_Bias - 1.96 * SE,
    CI_Upper = Mean_Bias + 1.96 * SE
  )

p_ci <- ggplot(df_ci, aes(x = Model, y = Mean_Bias, color = Model)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 1) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  facet_grid(Parameter ~ N_Label, scales = "free_y") +
  scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.8) +
  theme_n() +
  labs(title = "",
       y = expression("Mean Bias" %+-%"95% CI"), x = "Calibration Algorithm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

ggsave(
  here::here("outputs", "monte_carlo", "95pct_Confidence_Intervals_of_Estimation_Bias.pdf"), 
  plot = p_ci, 
  width = 14, 
  height = 8, 
  dpi = 300, 
  bg = "white"
)

# EXECUTION TIME WITH 95% CI
#-------------------------------------------------------------------------------
df_time <- all_metrics |> 
  dplyr::mutate(
    Model=dplyr::recode(Model, "F2017" = "BF")
  ) |>
  dplyr::group_by(Model, N_Obs) |> 
  dplyr::summarise(
    Mean_Time = mean(Time_sec, na.rm = TRUE),
    SD_Time = sd(Time_sec, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    CI_Lower = pmax(0.001, Mean_Time - 1.96 * (SD_Time / sqrt(N_sim))),
    CI_Upper = Mean_Time + 1.96 * (SD_Time / sqrt(N_sim))
  )

p_time <- ggplot(df_time, aes(x = N_Obs, y = Mean_Time, color = Model, fill = Model, group = Model)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.8) + 
  scale_fill_viridis_d(option = "mako", begin = 0.1, end = 0.8) + 
  scale_y_log10(labels = scales::comma) + 
  scale_x_continuous(breaks = unique(df_time$N_Obs)) +
  theme_n() +
  labs(title = "",
       x = "Number of Observations", y = "Execution Time (Seconds, Log Scale)")

ggsave(here::here("outputs", "monte_carlo", "Computational Scalability: Execution Time vs. Sample Siz.pdf"), plot = p_time, width = 9, height = 6, dpi = 300, bg = "white")

#RMSE WITH 95% CI RIBBONS
#-------------------------------------------------------------------------------
df_rmse <- all_metrics |> 
  dplyr::group_by(Model, N_Obs) |> 
  dplyr::summarise(
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    SD_RMSE = sd(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    CI_Lower = Mean_RMSE - 1.96 * (SD_RMSE / sqrt(N_sim)),
    CI_Upper = Mean_RMSE + 1.96 * (SD_RMSE / sqrt(N_sim)),
    Model=recode(
      Model, "F2017"="BF"
    )
  )

p_rmse <- ggplot(df_rmse, aes(x = N_Obs, y = Mean_RMSE, color = Model, fill = Model, group = Model)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
  scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.8) + 
  scale_fill_viridis_d(option = "mako", begin = 0.1, end = 0.8) + 
  scale_x_continuous(breaks = unique(df_rmse$N_Obs)) +
  theme_n() +
  labs(title = "",
       x = "Number of Observations", y = "Mean RMSE")

ggsave(here::here("outputs", "monte_carlo", "Model Accuracy: Root Mean Square Error (RMSE) vs. Sample Size.pdf"), plot = p_rmse, width = 9, height = 6, dpi = 300, bg = "white")

# PARETO FRONT
#-------------------------------------------------------------------------------
metrics_summary <- all_metrics |> 
  dplyr::mutate(
    Model = dplyr::recode(Model, "F2017" = "BF")
  ) |>
  dplyr::group_by(Model, N_Obs) |> 
  dplyr::summarise(
    Mean_RMSE        = mean(RMSE, na.rm = TRUE),
    Mean_Time_sec    = mean(Time_sec, na.rm = TRUE),
    Mean_Watson_U2   = mean(U2, na.rm = TRUE),
    Mode_Neyman_S    = as.numeric(names(sort(table(S_star), decreasing=TRUE)[1])), 
    Mean_TB2_Entropy = mean(TB2, na.rm = TRUE),
    .groups          = "drop"
  ) |> 
  dplyr::arrange(N_Obs, Mean_RMSE) 

df_tradeoff <- metrics_summary 

p_tradeoff <- ggplot(df_tradeoff, aes(x = Mean_Time_sec, y = Mean_RMSE)) +
  geom_point(aes(fill = Model), shape = 21, color = "black", 
             size = 5, stroke = 0.6, alpha = 0.85) +
  geom_text_repel(aes(label = Model, color = Model), 
                  size = 4.5, fontface = "bold", 
                  box.padding = 0.8, point.padding = 0.5,
                  segment.color = "gray50", 
                  show.legend = FALSE) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10() +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
  scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  ) + 
  labs(
    title = "",
    x = "Execution Time (Seconds, Log Scale)", 
    y = "Root Mean Square Error (Log Scale)"
  )
ggsave(here::here("outputs", "monte_carlo", "Fig_5_Comparative_Tradeoff.pdf"), plot = p_tradeoff, width = 9, height = 7, dpi = 300, bg = "white")

print(knitr::kable(metrics_summary, digits = 4, caption = "Table 1: Diagnostics for Comparative Review Article"))
readr::write_csv(metrics_summary, here::here("outputs", "monte_carlo", "Table_1_Comparative_Diagnostics.csv"))

# --------------------------------------------------------------------------------
# REAL DATA APPLICATION (SÃO PAULO REAL DATA) - USING TDS MODEL
# --------------------------------------------------------------------------------
data_path <- here::here("data", "sao_paulo_mgtwr_panel_dual_geom.rds")
censobr::set_censobr_cache_dir(path = here::here("data", "censobr_cache"))

sp_code <- 3550308 
projected_crs <- 31983 
panel_years <- c(2000, 2010, 2022) 

if (!file.exists(data_path)) {
  sf_train <- load_local_spatial_data(here::here("Data", "geoportal_estacao_trem_v2", "estacao_trem_v2.shp"), projected_crs) |> 
    dplyr::select(nm_estacao, geometry)
  sf_bus <- load_local_spatial_data(here::here("Data", "geoportal_terminal_onibus_v2", "terminal_onibus_v2.shp"), projected_crs) |> 
    dplyr::select(nm_termina, geometry)
  sf_slum <- load_local_spatial_data(here::here("Data", "habita2geosampa_habi_favela2geosampa", "habi_favela2geosampa.shp"), projected_crs) |> 
    dplyr::select(nome, geometry)
  sf_tenement <- load_local_spatial_data(here::here("Data", "habita2geosampa_habi_cortico2geosampa", "habi_cortico2geosampa.shp"), projected_crs) |> 
    dplyr::select(nome, geometry)
  sf_irregular_lot <- load_local_spatial_data(here::here("Data", "habita2geosampa_habi_loteamento2geosampa", "habi_loteamento2geosampa.shp"), projected_crs) |> 
    dplyr::select(nome, geometry)
  
  sp_districts <- sf::read_sf(here::here("Data", "Bairros_Distritos_CidadeSP-20260106T021204Z-1-001", "Bairros_Distritos_CidadeSP", "LAYER_DISTRITO", "DEINFO_DISTRITO.shp")) |> 
    sf::st_transform(projected_crs) |> 
    sf::st_make_valid()
  
  mgtwr_panel <- purrr::map_dfr(panel_years, ~process_panel_year(.x, sp_code, projected_crs, sf_train, sf_bus, sf_slum, sf_tenement, sf_irregular_lot))
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------==
  # INFLATION ADJUSTMENT (IPCA DEFLATOR)
  # Deflating 2000 and 2010 nominal income to 2022 Brazilian Reais (BRL) base.
  # ------------------------------------------------------------------------------------------------------------------------------------------------==
  mgtwr_panel <- mgtwr_panel |>
    dplyr::mutate(
      avg_income = dplyr::case_when(
        census_year == 2000 ~ avg_income * 3.784, 
        census_year == 2010 ~ avg_income * 1.955, 
        TRUE ~ avg_income # 2022 values remain unchanged (1.000)
      ),
      log_income = log(avg_income + 1) 
    )
  
  mgtwr_panel <- mgtwr_panel |> tidyr::drop_na(log_density, log_income) |> sf::st_make_valid()
  
  mgtwr_centroids <- suppressWarnings(sf::st_centroid(mgtwr_panel))
  mgtwr_with_dist <- sf::st_join(mgtwr_centroids, sp_districts["NOME_DIST"], join = sf::st_intersects)
  
  mgtwr_panel$NOME_DIST <- mgtwr_with_dist$NOME_DIST
  
  dist_medians_2010 <- mgtwr_panel |>
    sf::st_drop_geometry() |>
    dplyr::filter(census_year == 2010) |>
    dplyr::group_by(NOME_DIST) |>
    dplyr::summarise(
      dist_med_apt = median(prop_apt, na.rm = TRUE) * 1.05,   
      dist_med_sew = median(prop_sewage, na.rm = TRUE) * 1.02, 
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dist_med_apt = pmin(dist_med_apt, 100), 
      dist_med_sew = pmin(dist_med_sew, 100)  
    )
  
  global_apt_2010 <- min(median(mgtwr_panel$prop_apt[mgtwr_panel$census_year == 2010], na.rm = TRUE) * 1.05, 100)
  global_sew_2010 <- min(median(mgtwr_panel$prop_sewage[mgtwr_panel$census_year == 2010], na.rm = TRUE) * 1.02, 100)
  
  mgtwr_panel <- mgtwr_panel |>
    dplyr::left_join(dist_medians_2010, by = "NOME_DIST") |>
    dplyr::mutate(
      prop_apt = dplyr::case_when(
        census_year == 2022 & is.na(prop_apt) & !is.na(dist_med_apt) ~ dist_med_apt,
        census_year == 2022 & is.na(prop_apt) & is.na(dist_med_apt) ~ global_apt_2010,
        TRUE ~ prop_apt
      ),
      prop_sewage = dplyr::case_when(
        census_year == 2022 & is.na(prop_sewage) & !is.na(dist_med_sew) ~ dist_med_sew,
        census_year == 2022 & is.na(prop_sewage) & is.na(dist_med_sew) ~ global_sew_2010,
        TRUE ~ prop_sewage
      )
    ) |>
    dplyr::select(-dist_med_apt, -dist_med_sew)
  
  sp_geom_df <- tibble::as_tibble(sp_districts) |> dplyr::select(NOME_DIST, geom_district = geometry)
  mgtwr_panel_dual <- mgtwr_panel |> dplyr::left_join(sp_geom_df, by = "NOME_DIST")
  sf::st_geometry(mgtwr_panel_dual) <- "geom"
  mgtwr_panel_dual <- sf::st_cast(mgtwr_panel_dual, "MULTIPOLYGON")
  
  saveRDS(mgtwr_panel_dual, data_path)
} else {
  mgtwr_panel_dual <- readRDS(data_path)
}

mgtwr_panel <- mgtwr_panel_dual |> dplyr::mutate(census_year = as.factor(census_year))
df_tabular <- mgtwr_panel |> sf::st_drop_geometry() |> dplyr::select(-dplyr::any_of(c("geom_district")))

districts_clean <- mgtwr_panel |> 
  sf::st_drop_geometry() |>
  dplyr::select(NOME_DIST, geom_district) |>
  dplyr::distinct(NOME_DIST, .keep_all = TRUE) |>
  sf::st_as_sf(sf_column_name = "geom_district")

sf::st_crs(districts_clean) <- sf::st_crs(mgtwr_panel)

# --------------------------------------------------------------------------------
# PART 4: ESDA (Exploratory Spatial Data Analysis)
# --------------------------------------------------------------------------------
# DESCRIPTIVE STATISTICS
desc_table <- df_tabular |>
  dplyr::select(census_year, density_sqkm, log_density, avg_income, log_income, prop_apt,
                prop_sewage, persons_per_hh, dist_train_km, dist_bus_km, 
                pct_area_slum, pct_area_tenement, pct_area_irregular_lot) |>
  gtsummary::tbl_summary(
    by         = census_year,
    statistic  = list(gtsummary::all_continuous() ~ "{mean} ({sd})"),
    digits     = gtsummary::all_continuous() ~ 2,
    label      = list(
      density_sqkm           ~ "Population Density (pop/sq.km)",
      log_density            ~ "Log Population Density (Y)",
      avg_income             ~ "Average Nominal Income (BRL)",
      log_income             ~ "Log Average Income",
      prop_apt               ~ "Households in Apartments (%)",
      prop_sewage            ~ "Adequate Sewage Coverage (%)",
      persons_per_hh         ~ "Household Size (Persons/HH)",
      dist_train_km          ~ "Distance to Train/Subway (km)",
      dist_bus_km            ~ "Distance to Bus Terminal (km)",
      pct_area_slum          ~ "Area overlapping Slums (%)",
      pct_area_tenement      ~ "Area overlapping Tenements (%)",
      pct_area_irregular_lot ~ "Area overlapping Irregular Lots (%)"
    )
  ) |>
  gtsummary::add_overall() |>
  gtsummary::modify_header(label ~ "**Variable**") |>
  gtsummary::modify_caption("**Table 2. Empirical Descriptive Statistics**")

desc_table |> gtsummary::as_gt() |> gt::gtsave(here::here("outputs", "real", "Table2_Descriptive_Statistics.docx"))

# CORRELATION MATRIX
cor_vars <- df_tabular |>
  dplyr::select(log_density, log_income, prop_apt, prop_sewage, persons_per_hh, 
                dist_train_km, dist_bus_km, pct_area_slum, pct_area_tenement, pct_area_irregular_lot) |>
  dplyr::rename(`Log Density` = log_density, `Log Income` = log_income,
                `Apts (%)` = prop_apt, `Sewage (%)` = prop_sewage,
                `HH Size` = persons_per_hh, `Dist. Train` = dist_train_km,
                `Dist. Bus` = dist_bus_km, `Slum (%)` = pct_area_slum,
                `Tenement (%)` = pct_area_tenement, `Irregular Lot (%)` = pct_area_irregular_lot)

cor_matrix <- cor(cor_vars, use = "complete.obs", method = "spearman")
png(here::here("outputs", "real", "Fig6_Correlation_Matrix.png"), width = 2800, height = 2800, res = 300)
corrplot::corrplot(cor_matrix, method = "color", type = "lower", tl.col = "black", tl.srt = 45,
                   addCoef.col = "black", number.cex = 0.7,
                   col = colorRampPalette(c("#440154FF", "white", "#FDE725FF"))(200),
                   title = "Spearman Correlation Matrix", mar = c(0, 0, 1, 0))
dev.off()

# DISTRIBUTIONS
ridge_density <- ggplot(df_tabular, aes(x = log_density, y = census_year, fill = census_year)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", show.legend = FALSE) +
  scale_fill_viridis_d(option = "mako") + theme_n() + 
  labs(title = "A. Population Density Evolution", x = "Log(Density)", y = "Census Year")

ridge_income <- ggplot(df_tabular, aes(x = log_income, y = census_year, fill = census_year)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", show.legend = FALSE) +
  scale_fill_viridis_d(option = "mako") + theme_n() + 
  labs(title = "B. Income Distribution Evolution", x = "Log(Income)", y = "")

fig_distribution <- ridge_density + ridge_income
ggsave(here::here("outputs", "real", "Fig7_Distributions.png"), plot = fig_distribution, width = 12, height = 5, dpi = 300, bg = "white")

# DENSITY MAP
district_density_panel <- df_tabular |>
  dplyr::group_by(census_year, NOME_DIST) |>
  dplyr::summarise(
    log_density = mean(log_density, na.rm = TRUE), 
    .groups = "drop"
  ) |>
  dplyr::left_join(districts_clean, by = "NOME_DIST") |>
  sf::st_as_sf(sf_column_name = "geom_district")

# --------------------------------------------------------------------------------
# DENSITY MAP
# --------------------------------------------------------------------------------

map_density <- ggplot2::ggplot(district_density_panel) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = log_density), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::scale_fill_fermenter(
    palette = "YlOrRd", 
    direction = 1, 
    name = "Log(Density)",
    breaks = scales::pretty_breaks(n = 7), 
    guide = ggplot2::guide_coloursteps(
      barwidth = ggplot2::unit(12, "lines"),
      barheight = ggplot2::unit(0.5, "lines"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + 
  map_theme_addons +
  ggplot2::labs(
    title = ""
  ) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 13, hjust = 0.5, color = "grey20", margin = ggplot2::margin(b = 15)),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  here::here("outputs", "real", "Fig8_Map_Density_Publication.pdf"), 
  plot = map_density, width = 15, height = 7, device = "pdf"
)

ggplot2::ggsave(here::here("outputs", "real", "Fig8_Map_Density.png"), plot = map_density, width = 14, height = 6, dpi = 300, bg = "white")

# --------------------------------------------------------------------------------
# LISA MAPS & SPATIAL AUTOCORRELATION
# --------------------------------------------------------------------------------

lisa_panel <- map(c("2000", "2010", "2022"), \(yr) calculate_lisa_by_year(yr, sf_data = mgtwr_panel)) |> 
  list_rbind() |> 
  st_as_sf() 

st_geometry(lisa_panel) <- st_geometry(lisa_panel) 
st_crs(lisa_panel) <- st_crs(districts_clean) 

lisa_colors <- c(
  "High-High (HH)"  = "#e31a1c", 
  "Low-Low (LL)"    = "#1f78b4", 
  "High-Low (HL)"   = "#fb9a99", 
  "Low-High (LH)"   = "#a6cee3", 
  "Not Significant" = "#eeeeee" 
)

lisa_map <- ggplot2::ggplot(lisa_panel) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = cluster), color = "black", linewidth = 0.1, alpha = 0.9) +
  
  ggplot2::scale_fill_manual(
    values = lisa_colors, 
    name = "LISA Cluster Type",
    guide = ggplot2::guide_legend(
      nrow = 1, 
      title.position = "top", 
      title.hjust = 0.5,
      label.position = "bottom",
      keywidth = 3
    )
  ) +
  
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + 
  map_theme_addons +
  ggplot2::labs(
    title = ""
  ) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 13, hjust = 0.5, color = "grey20", margin = ggplot2::margin(b = 15)),
    legend.position = "bottom",
    strip.text    = ggplot2::element_text(size = 12, face = "bold"),
    panel.border  = ggplot2::element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

ggplot2::ggsave(
  here::here("outputs", "real", "Fig4_LISA_Spatiotemporal_Map.pdf"), 
  plot = lisa_map, width = 15, height = 7, device = "pdf"
)

ggplot2::ggsave(
  here::here("outputs", "real", "Fig4_LISA_Spatiotemporal_Map.png"), 
  plot = lisa_map, width = 14, height = 6, dpi = 300, bg = "white"
)

stats_table <- lisa_panel |> 
  sf::st_drop_geometry() |> 
  dplyr::group_by(census_year) |> 
  dplyr::summarise(
    Moran_I_Numeric = dplyr::first(moran_global_i),
    P_Value_Numeric = dplyr::first(moran_p_value),
    .groups = "drop"
  ) |> 
  dplyr::mutate(
    `Global Moran's I` = round(Moran_I_Numeric, 4),
    `P-Value`          = format.pval(P_Value_Numeric, digits = 3, eps = 0.001),
    Interpretation     = dplyr::case_when(
      P_Value_Numeric < 0.05 & Moran_I_Numeric > 0 ~ "Significant POSITIVE Autocorrelation",
      P_Value_Numeric < 0.05 & Moran_I_Numeric < 0 ~ "Significant NEGATIVE Autocorrelation",
      TRUE                                         ~ "Random (Not Significant)"
    )
  ) |> 
  dplyr::select(census_year, `Global Moran's I`, `P-Value`, Interpretation) |>
  dplyr::rename(`Census Year` = census_year)

knitr::kable(stats_table, 
             align   = c("c", "c", "c", "l"),
             caption = "Table 3. Evolution of Global Spatial Autocorrelation (Moran's I) for Population Density") |>
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) |>
  print()

# --------------------------------------------------------------------------------
# TDS-MGTWR
# --------------------------------------------------------------------------------
model_data_real <- mgtwr_panel |>
  tidyr::drop_na(log_income, log_density, prop_apt, prop_sewage,
                 dist_train_km, dist_bus_km, pct_area_slum, persons_per_hh, 
                 pct_area_irregular_lot, pct_area_tenement) |>
  sf::st_as_sf() 

sf::st_geometry(model_data_real) <- "geom"

sd_dict <- c(
  "Z_income"     = sd(model_data_real$log_income, na.rm = TRUE),
  "Z_apt"        = sd(log(model_data_real$prop_apt + 1), na.rm = TRUE),
  "Z_sewage"     = sd(model_data_real$prop_sewage, na.rm = TRUE),
  "Z_train"      = sd(model_data_real$dist_train_km, na.rm = TRUE),
  "Z_bus"        = sd(model_data_real$dist_bus_km, na.rm = TRUE),
  "Z_slum"       = sd(log(model_data_real$pct_area_slum + 1), na.rm = TRUE),
  "Z_tenement"   = sd(log(model_data_real$pct_area_tenement + 1), na.rm = TRUE),
  "Z_irregular"  = sd(log(model_data_real$pct_area_irregular_lot + 1), na.rm = TRUE),
  "Z_persons_hh" = sd(model_data_real$persons_per_hh, na.rm = TRUE)
)
sd_y <- sd(model_data_real$log_density, na.rm = TRUE)

# Standardize variables (Z-scores)
model_data_real <- model_data_real |>
  dplyr::mutate(
    Z_density    = scale(log_density)[, 1],
    Z_income     = scale(log_income)[, 1],
    Z_apt        = scale(log(prop_apt + 1))[, 1],
    Z_sewage     = scale(prop_sewage)[, 1],
    Z_train      = scale(dist_train_km)[, 1],
    Z_bus        = scale(dist_bus_km)[, 1],
    Z_slum       = scale(log(pct_area_slum + 1))[, 1],
    Z_tenement   = scale(log(pct_area_tenement + 1))[, 1],
    Z_irregular  = scale(log(pct_area_irregular_lot + 1))[, 1],
    Z_persons_hh = scale(persons_per_hh)[, 1]
  )

formula_real <- Z_density ~ Z_income + Z_apt + Z_sewage + Z_train + Z_bus + Z_slum + Z_tenement + Z_irregular + Z_persons_hh

if (requireNamespace("car", quietly = TRUE)) {
  global_ols <- lm(formula_real, data = model_data_real)
  
  vif_vals <- car::vif(global_ols)
  vif_df <- data.frame(Variable = names(vif_vals), VIF = as.numeric(vif_vals))
  readr::write_csv(vif_df, here::here("outputs", "real", "Table_VIF_Diagnostics.csv"))
  
  if (requireNamespace("spdep", quietly = TRUE)) {
    coords_ols <- sf::st_coordinates(suppressWarnings(sf::st_centroid(model_data_real)))
    knn_ols <- spdep::knearneigh(coords_ols, k = 5)
    nb_ols <- spdep::knn2nb(knn_ols)
    listw_ols <- spdep::nb2listw(nb_ols, style = "W", zero.policy = TRUE)
    
    ols_resid_moran <- spdep::moran.test(resid(global_ols), listw_ols, zero.policy = TRUE)
    
    ols_moran_df <- data.frame(
      Diagnostic = "OLS Residual Moran's I",
      Statistic = ols_resid_moran$estimate[1],
      P_Value = ols_resid_moran$p.value
    )
    readr::write_csv(ols_moran_df, here::here("outputs", "real", "Table_OLS_Residual_Moran.csv"))
  }
}

# SAMPLING
set.seed(2026)
model_data_final <- model_data_real %>%
  dplyr::group_by(census_year, NOME_DIST) %>%
  dplyr::slice_sample(n = 250, replace = FALSE) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(census_year) %>%
  dplyr::slice_sample(n = 15000) %>% 
  dplyr::ungroup()

# GLOBAL OLS VALIDATION
ols_validation <- lm(formula_real, data = model_data_final)
ols_estimates <- data.frame(
  Variable = names(coef(ols_validation)),
  Global_OLS = as.numeric(coef(ols_validation))
) %>%
  dplyr::mutate(Variable = ifelse(Variable == "(Intercept)", "Intercept", Variable))

# TDS-MGTWR MODEL EXECUTION
res_real <- evaluate_and_plot_mgtwr(
  model_data      = model_data_final,
  districts_clean = districts_clean,
  prefix_name     = "SaoPaulo_Empirical",
  formula         = formula_real,
  models_to_run   = "TDS",  
  out_dir         = here::here("outputs", "real")
)

readr::write_csv(res_real$metrics, here::here("outputs", "real", "Table_MGTWR_Empirical_Metrics.csv"))

if (!is.null(res_real$mod_tds$optimal_bws)) {
  x_vars <- attr(terms(formula_real), "term.labels")
  bws_mat <- res_real$mod_tds$optimal_bws
  
  bws_df <- data.frame(
    Variable = c("Intercept", x_vars),
    Spatial_Bandwidth = bws_mat[, 1],
    Temporal_Bandwidth = bws_mat[, 2]
  )
  readr::write_csv(bws_df, here::here("outputs", "real", "Table_TDS_Bandwidths.csv"))
}

model_data_real_out <- res_real$model_data |> sf::st_as_sf()
sf::st_geometry(model_data_real_out) <- "geom"

model_data_real_out$yhat <- res_real$mod_tds$yhat
model_data_real_out$residual <- res_real$mod_tds$residual

x_vars <- attr(terms(formula_real), "term.labels")
model_data_real_out[["Beta_Intercept_Original"]] <- model_data_real_out[["Beta_Intercept_tds"]]

for (var in x_vars) {
  beta_col <- paste0("Beta_", var, "_tds")
  orig_col <- paste0("Beta_", var, "_Original")
  sd_x <- sd_dict[[var]]
  model_data_real_out[[orig_col]] <- model_data_real_out[[beta_col]] * (sd_y / sd_x)
}

local_coefs_summary <- res_real$model_data %>%
  sf::st_drop_geometry() %>%
  dplyr::select(starts_with("Beta_") & ends_with("_tds")) %>%
  tidyr::pivot_longer(everything(), names_to = "Variable", values_to = "Coefficient") %>%
  dplyr::mutate(Variable = stringr::str_replace_all(Variable, c("Beta_" = "", "_tds" = ""))) %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(
    Mean     = mean(Coefficient, na.rm = TRUE),
    SD       = sd(Coefficient, na.rm = TRUE),
    Min      = min(Coefficient, na.rm = TRUE),
    P25      = quantile(Coefficient, 0.25, na.rm = TRUE),
    Median   = median(Coefficient, na.rm = TRUE),
    P75      = quantile(Coefficient, 0.75, na.rm = TRUE),
    Max      = max(Coefficient, na.rm = TRUE),
    IQ_Range = P75 - P25,
    .groups  = "drop"
  )

final_heterogeneity_table <- local_coefs_summary %>%
  dplyr::left_join(ols_estimates, by = "Variable") %>%
  dplyr::select(Variable, Global_OLS, Mean, Median, SD, Min, Max, IQ_Range) %>%
  dplyr::mutate(Spatial_Range = abs(Max - Min)) %>%
  dplyr::arrange(desc(Spatial_Range))

print("--- SPATIOTEMPORAL HETEROGENEITY SUMMARY ---")
print(final_heterogeneity_table)

readr::write_csv(final_heterogeneity_table, here::here("outputs", "real", "Table_Local_Beta_Distribution.csv"))

global_ols <- lm(formula_real, data = model_data_real)

ols_coefs <- data.frame(Variable = names(coef(global_ols)), Global_Estimate = as.numeric(coef(global_ols))) |> 
  dplyr::mutate(Variable = stringr::str_replace(Variable, "\\(Intercept\\)", "Intercept")) |>
  dplyr::mutate(Variable = paste0("Beta_", Variable, "_tds"))

df_local_coefs <- model_data_real_out |> 
  sf::st_drop_geometry() |> 
  dplyr::select(starts_with("Beta_") & ends_with("_tds")) |>
  tidyr::pivot_longer(everything(), names_to = "Variable", values_to = "Local_Estimate") |>
  dplyr::left_join(ols_coefs, by = "Variable") |>
  dplyr::mutate(
    Variable_Clean = stringr::str_replace_all(Variable, c("Beta_Z_" = "", "_tds" = "", "Beta_" = "")),
    Variable_Clean = stringr::str_to_title(Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Persons_hh", "Household Size", Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Dist_train", "Dist. Train", Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Dist_bus", "Dist. Bus", Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Slum", "Slums Area", Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Tenement", "Tenements Area", Variable_Clean),
    Variable_Clean = ifelse(Variable_Clean == "Irregular", "Irregular Lots Area", Variable_Clean)
  )

df_forest_summary <- df_local_coefs |>
  dplyr::group_by(Variable_Clean, Global_Estimate) |>
  dplyr::summarise(
    Spatial_Median = median(Local_Estimate, na.rm = TRUE),
    Spatial_P025 = quantile(Local_Estimate, 0.025, na.rm = TRUE),
    Spatial_P975 = quantile(Local_Estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_forest <- ggplot(df_local_coefs, aes(x = Local_Estimate, y = reorder(Variable_Clean, Global_Estimate))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) +
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.2, point_colour = NA, fill = "slategray", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", alpha = 0.8) +
  geom_point(data = df_forest_summary, aes(x = Global_Estimate, y = Variable_Clean, fill = "Global OLS"), shape = 23, size = 4, color = "darkred", fill = "red") +
  theme_n() +
  labs(title = "Spatial Heterogeneity of Determinants", 
       x = "Standardized Coefficient Value", y = "Predictor Variable")

ggsave(here::here("outputs", "real", "Fig9_Spatial_Heterogeneity_Forest.png"), plot = p_forest, width = 10, height = 7, dpi = 300, bg = "white")


full_map_data <- model_data_real |> 
  sf::st_as_sf() |> 
  sf::st_set_geometry("geom")

full_map_data <- full_map_data[, !(names(full_map_data) %in% "geom_district")]

sampled_centroids <- model_data_real_out |> 
  sf::st_as_sf() |> 
  sf::st_set_geometry("geom") |> 
  sf::st_centroid()

interpolated_years <- list()

for (yr in unique(full_map_data$census_year)) {
  full_yr <- full_map_data |> dplyr::filter(census_year == yr)
  samp_yr <- sampled_centroids |> dplyr::filter(census_year == yr)
  full_yr <- sf::st_as_sf(full_yr)
  
  joined_yr <- full_yr |>
    sf::st_join(
      samp_yr |> dplyr::select(starts_with("Beta_"), yhat, residual), 
      join = sf::st_nearest_feature
    )
  
  interpolated_years[[as.character(yr)]] <- joined_yr
}

model_data_plot_full <- dplyr::bind_rows(interpolated_years) |> sf::st_as_sf()

model_data_plot_full <- model_data_plot_full |> 
  dplyr::mutate(
    Beta_Income     = Beta_Z_income_Original,
    Beta_Apt        = Beta_Z_apt_Original,
    Beta_Sewage     = Beta_Z_sewage_Original,
    Beta_Train      = Beta_Z_train_Original,
    Beta_Slum       = Beta_Z_slum_Original,
    Beta_Bus        = Beta_Z_bus_Original,
    Beta_Tenement   = Beta_Z_tenement_Original,
    Beta_Irregular  = Beta_Z_irregular_Original,
    Beta_Persons_hh = Beta_Z_persons_hh_Original,
    Beta_Intercept  = Beta_Intercept_Original
  )

path_out <- here::here("outputs", "real")

# ------------------------------------------------------------------------------
# CALCULATING NATURAL BREAKS (fisher) TO AVOID DEFAULT DISTORTION
# ------------------------------------------------------------------------------
get_fisher <- function(variable, n_classes = 6) {
  valid_data <- na.omit(variable)
  if(length(unique(valid_data)) < n_classes) return(scales::pretty_breaks(n = n_classes)(valid_data)) 
  suppressWarnings(classInt::classIntervals(valid_data, n = n_classes, style = "fisher")$brks)
}

breaks_yhat           <- get_fisher(model_data_plot_full$yhat)
breaks_raw_income     <- get_fisher(model_data_plot_full$Beta_Income)
breaks_raw_apt        <- get_fisher(model_data_plot_full$Beta_Apt)
breaks_raw_sewage     <- get_fisher(model_data_plot_full$Beta_Sewage)
breaks_raw_train      <- get_fisher(model_data_plot_full$Beta_Train)
breaks_raw_slum       <- get_fisher(model_data_plot_full$Beta_Slum)
breaks_raw_bus        <- get_fisher(model_data_plot_full$Beta_Bus)
breaks_raw_tenement   <- get_fisher(model_data_plot_full$Beta_Tenement)
breaks_raw_irregular  <- get_fisher(model_data_plot_full$Beta_Irregular)
breaks_raw_persons_hh <- get_fisher(model_data_plot_full$Beta_Persons_hh)

fermenter_guide <- ggplot2::guide_colorsteps(
  barwidth = ggplot2::unit(15, "lines"), barheight = ggplot2::unit(0.5, "lines"),
  title.position = "top", title.hjust = 0.5, show.limits = TRUE
)

# ------------------------------------------------------------------------------
# ALL MAPS
# ------------------------------------------------------------------------------

# Predicted Log Density
p_yhat <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = yhat), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "YlOrRd", direction = 1, name = expression(hat(Y)), breaks = breaks_yhat, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Income
p_raw_income <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Income), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Income"]), breaks = breaks_raw_income, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Verticalization (Apts)
p_raw_apt <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Apt), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Apartments"]), breaks = breaks_raw_apt, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Sewage Coverage
p_raw_sewage <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Sewage), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Sewage"]), breaks = breaks_raw_sewage, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Mobility (Train)
p_raw_train <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Train), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Train"]), breaks = breaks_raw_train, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Informality (Slums)
p_raw_slum <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Slum), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = 1, name = expression(hat(beta)["Slum"]), breaks = breaks_raw_slum, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Mobility (Bus)
p_raw_bus <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Bus), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Bus"]), breaks = breaks_raw_bus, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Tenements
p_raw_tenement <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Tenement), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Tenement"]), breaks = breaks_raw_tenement, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Irregular Lots
p_raw_irregular <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Irregular), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Irregular"]), breaks = breaks_raw_irregular, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Local Effect: Household Size
p_raw_persons_hh <- ggplot2::ggplot(model_data_plot_full) +
  ggspatial::annotation_map_tile(type = "cartolight", progress = "none", alpha = 0.5) +
  ggplot2::geom_sf(ggplot2::aes(fill = Beta_Persons_hh), color = "white", linewidth = 0.1, alpha = 0.9) +
  ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_fermenter(palette = "RdBu", direction = -1, name = expression(hat(beta)["Household"]), breaks = breaks_raw_persons_hh, guide = fermenter_guide, labels = scales::label_number(accuracy = 0.01)) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.3) +
  ggspatial::annotation_north_arrow(location = "br", style = ggspatial::north_arrow_minimal(), pad_y=unit(0.55, "cm")) +
  ggplot2::facet_wrap(~ census_year, ncol = 3) +
  theme_n() + map_theme_addons +
  ggplot2::labs(title = "", x = NULL, y = NULL) +
  ggplot2::theme(legend.position = "bottom")

# Export all maps
ggplot2::ggsave(file.path(path_out, "Fig11_A_Yhat.png"), plot = p_yhat, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_B_Income_Raw.png"), plot = p_raw_income, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_C_Apt_Raw.png"), plot = p_raw_apt, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_D_Sewage_Raw.png"), plot = p_raw_sewage, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_E_Train_Raw.png"), plot = p_raw_train, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_F_Slum_Raw.png"), plot = p_raw_slum, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_G_Bus_Raw.png"), plot = p_raw_bus, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_H_Tenement_Raw.png"), plot = p_raw_tenement, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_I_Irregular_Raw.png"), plot = p_raw_irregular, width = 10, height = 7, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(path_out, "Fig11_J_Persons_hh_Raw.png"), plot = p_raw_persons_hh, width = 10, height = 7, dpi = 300, bg = "white")
