#-------------------------------------------------------------------------------
# SPATIO-TEMPORAL REGRESSION UTILITIES & ALGORITHMS
# Focus: Comparative Analysis of Existing MGTWR Calibration Methods
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# SPATIAL DATA UTILITIES
#-------------------------------------------------------------------------------

ensure_crs <- function(sf_obj, crs) {
  if (is.null(sf::st_crs(sf_obj))) sf_obj <- sf::st_set_crs(sf_obj, crs)
  if (sf::st_crs(sf_obj)$epsg != sf::st_crs(crs)$epsg) sf_obj <- sf::st_transform(sf_obj, crs)
  sf::st_make_valid(sf_obj)
}

safe_centroids <- function(sf_obj, method = c("point_on_surface", "centroid")) {
  method <- match.arg(method)
  if (method == "point_on_surface") {
    pts <- suppressWarnings(sf::st_point_on_surface(sf_obj))
  } else {
    pts <- suppressWarnings(sf::st_centroid(sf_obj))
  }
  return(sf::st_coordinates(pts))
}

load_local_spatial_data <- function(path, crs) {
  if (!file.exists(path)) stop(paste("File not found:", path))
  sf <- sf::read_sf(path)
  return(ensure_crs(sf, crs))
}

load_all_spatial_layers <- function(paths_named_list, crs) {
  res <- list()
  for (nm in names(paths_named_list)) {
    res[[nm]] <- tryCatch({
      load_local_spatial_data(paths_named_list[[nm]], crs)
    }, error = function(e) {
      message(sprintf("[WARN] Could not load %s: %s", nm, e$message))
      NULL
    })
  }
  return(res)
}

spatial_impute <- function(val_vector, nb_obj) {
  imputed <- val_vector
  na_idx  <- which(is.na(val_vector) | val_vector == 0)
  for (i in na_idx) {
    neighbors  <- nb_obj[[i]]
    neighbors  <- neighbors[neighbors != 0]
    if (length(neighbors) > 0) {
      neighbor_vals <- val_vector[neighbors]
      valid_vals    <- neighbor_vals[!is.na(neighbor_vals) & neighbor_vals > 0]
      if (length(valid_vals) > 0) imputed[i] <- mean(valid_vals)
    }
  }
  return(imputed)
}

calc_informality_overlap <- function(sf_tracts, sf_informal) {
  sf_informal_u <- suppressWarnings(sf::st_union(sf::st_make_valid(sf_informal)))
  isect <- suppressWarnings(sf::st_intersection(sf::st_make_valid(sf_tracts), sf_informal_u))
  
  if (nrow(isect) == 0) return(rep(0, nrow(sf_tracts)))
  
  isect$area_inf <- as.numeric(sf::st_area(isect)) / 1e6
  df_inf <- isect |> sf::st_drop_geometry() |>
    dplyr::group_by(code_tract) |>
    dplyr::summarise(area_inf = sum(area_inf), .groups = "drop")
  
  res <- sf_tracts |> sf::st_drop_geometry() |>
    dplyr::left_join(df_inf, by = "code_tract") |>
    dplyr::mutate(pct = (dplyr::coalesce(area_inf, 0) / (as.numeric(sf::st_area(sf_tracts)) / 1e6)) * 100)
  return(res$pct)
}

calculate_lisa_analysis <- function(sf_data, var_name) {
  nb <- spdep::poly2nb(sf_data, queen = TRUE)
  no_nb <- which(spdep::card(nb) == 0)
  
  if (length(no_nb) > 0) {
    coords <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(sf_data))))
    for (i in no_nb) {
      nb[[i]] <- as.integer(order(sf::st_distance(suppressWarnings(sf::st_centroid(sf_data[i, ])), suppressWarnings(sf::st_centroid(sf_data))))[2])
    }
  }
  
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  lisa <- spdep::localmoran(sf_data[[var_name]], listw)
  
  scaled_v <- as.numeric(scale(sf_data[[var_name]]))
  lag_v <- as.numeric(spdep::lag.listw(listw, scaled_v))
  
  sf_data$cluster <- dplyr::case_when(
    scaled_v > 0 & lag_v > 0 & lisa[, 5] < 0.05 ~ "High-High (HH)",
    scaled_v < 0 & lag_v < 0 & lisa[, 5] < 0.05 ~ "Low-Low (LL)",
    scaled_v > 0 & lag_v < 0 & lisa[, 5] < 0.05 ~ "High-Low (HL)",
    scaled_v < 0 & lag_v > 0 & lisa[, 5] < 0.05 ~ "Low-High (LH)",
    TRUE ~ "Not Significant"
  )
  
  sf_data$moran_global_i <- spdep::moran.test(sf_data[[var_name]], listw = listw, randomisation = TRUE, zero.policy = TRUE)$estimate[1]
  sf_data$moran_p_value  <- spdep::moran.test(sf_data[[var_name]], listw = listw, randomisation = TRUE, zero.policy = TRUE)$p.value
  
  sf_data <- sf_data |> dplyr::mutate(cluster = factor(cluster, levels = c("High-High (HH)", "Low-Low (LL)", "High-Low (HL)", "Low-High (LH)", "Not Significant")))
  return(sf_data)
}

#-------------------------------------------------------------------------------
# CENSUS EXTRACTION & HARMONIZATION
#-------------------------------------------------------------------------------

safe_census_download <- function(year, dataset_options, sp_code) {
  for (ds in dataset_options) {
    data <- tryCatch({
      res <- censobr::read_tracts(year = year, dataset = ds, showProgress = FALSE) |> 
        dplyr::collect()
      
      names(res) <- toupper(names(res))
      
      muni_col <- intersect(c("CODE_MUNI", "COD_MUNI", "COD_MUNICIPIO", "MUNICIPIO"), names(res))[1]
      tract_col <- intersect(c("CODE_TRACT", "COD_SETOR", "SETOR"), names(res))[1]
      
      if (!is.na(tract_col)) {
        res[[tract_col]] <- as.character(res[[tract_col]])
        
        if (!is.na(muni_col)) {
          res <- res[res[[muni_col]] == sp_code | res[[muni_col]] == as.character(sp_code), ]
        } else if (year == 2010) {
          res <- res[substr(res[[tract_col]], 1, nchar(as.character(sp_code))) == as.character(sp_code), ]
        }
        
        res <- res[!duplicated(res[[tract_col]]), ]
        names(res)[names(res) == tract_col] <- "CODE_TRACT"
        res
      } else {
        NULL
      }
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(data) && nrow(data) > 0) {
      message(sprintf("    [OK] Loaded universe dataset: %s", ds))
      return(data)
    }
  }
  message(sprintf(" > Warning: Dataset not found or empty for %d (Tried: %s)", year, paste(dataset_options, collapse=", ")))
  return(NULL)
}

extract_column_multi <- function(df_main, df_alt, new_col_name, possible_names, df_ref_hh = NULL, min_prop = 0, max_prop = 1) {
  check_and_return <- function(df, found_col, source_name) {
    val_vec <- as.numeric(df[[found_col]])
    if (!is.null(df_ref_hh)) {
      ref_vec <- as.numeric(df_ref_hh)
      valid_idx <- which(ref_vec > 0 & !is.na(val_vec))
      if (length(valid_idx) > 0) {
        prop <- mean(val_vec[valid_idx] / ref_vec[valid_idx], na.rm = TRUE)
        if (prop < min_prop || prop > max_prop) {
          return(NULL)
        }
      }
    }
    df_selected <- df |> dplyr::select(CODE_TRACT, dplyr::all_of(found_col))
    names(df_selected)[2] <- new_col_name
    return(df_selected)
  }
  
  if (!is.null(df_main)) {
    for (col in intersect(possible_names, names(df_main))) {
      res <- check_and_return(df_main, col, "Main DB")
      if (!is.null(res)) return(res)
    }
  }
  if (!is.null(df_alt)) {
    for (col in intersect(possible_names, names(df_alt))) {
      res <- check_and_return(df_alt, col, "Alt DB")
      if (!is.null(res)) return(res)
    }
  }
  if (is.null(df_main)) return(tibble::tibble(CODE_TRACT = character(0), !!new_col_name := numeric(0)))
  return(df_main |> dplyr::select(CODE_TRACT) |> dplyr::mutate(!!new_col_name := NA_real_))
}

process_panel_year <- function(year, sp_code, projected_crs, sf_train, sf_bus, sf_slum, sf_tenement, sf_irregular_lot) {
  message(sprintf("\nStarting Processing for Census Year %d...", year))
  
  sf_tracts <- geobr::read_census_tract(code_tract = sp_code, year = year, showProgress = FALSE) |>
    sf::st_transform(projected_crs) |> sf::st_make_valid() |>
    dplyr::mutate(code_tract = as.character(code_tract)) |>
    dplyr::select(code_tract, dplyr::any_of(c("geom", "geometry")))
  
  # Download essential census tables
  if (year == 2000) {
    df_basic     <- safe_census_download(year, c("Basico", "basico"), sp_code)
    df_household <- safe_census_download(year, c("Domicilio", "domicilio", "Domicilio01"), sp_code)
    df_persons   <- safe_census_download(year, c("Pessoa", "pessoa", "Pessoas", "pessoas", "Pessoa03"), sp_code)
    
    p_pop <- extract_column_multi(df_basic, df_persons, "total_pop", c("VAR12", "V002"))
    p_hh  <- extract_column_multi(df_basic, df_household, "total_hh",  c("VAR01", "V0003", "V001"))
    p_inc <- extract_column_multi(df_basic, df_persons, "avg_income", c("VAR06", "V005"))
    d_apt <- extract_column_multi(df_household, NULL, "hh_apt",    c("V0007", "V007"),    df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
    d_sew <- extract_column_multi(df_household, NULL, "hh_sewage", c("V0030", "V030"),    df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
    
  } else if (year == 2010) {
    df_basic     <- safe_census_download(year, c("Basico", "basico"), sp_code)
    df_household <- safe_census_download(year, c("Domicilio01", "domicilio01", "Domicilio", "domicilios"), sp_code)
    
    p_pop <- extract_column_multi(df_basic, NULL, "total_pop", c("V002", "BASICO_V0002"))
    p_hh  <- extract_column_multi(df_basic, NULL, "total_hh",  c("V001", "BASICO_V0001"))
    p_inc <- extract_column_multi(df_basic, NULL, "avg_income", c("V005", "BASICO_V0005"))
    d_apt <- extract_column_multi(df_household, NULL, "hh_apt",    c("V005", "DOMICILIO01_V005", "V0005"), df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
    d_sew <- extract_column_multi(df_household, NULL, "hh_sewage", c("V017", "DOMICILIO01_V017", "V0017"), df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
    
  } else if (year == 2022) {
    df_basic     <- safe_census_download(year, c("Basico", "basico"), sp_code)
    df_household <- safe_census_download(year, c("Domicilio", "domicilios"), sp_code)
    df_persons   <- safe_census_download(year, c("Pessoas", "pessoas", "Pessoa", "pessoa"), sp_code)
    df_income    <- safe_census_download(year, c("ResponsavelRenda", "responsavelrenda", "Renda", "renda"), sp_code)
    
    p_pop <- extract_column_multi(df_basic, df_persons, "total_pop", c("V0001", "V001"))
    p_hh  <- extract_column_multi(df_basic, df_household, "total_hh",  c("V0002", "V002"))
    p_inc <- extract_column_multi(df_income, df_basic, "avg_income", c("V06004", "V0005"))
    d_apt <- extract_column_multi(df_household, NULL, "hh_apt",    c("V0008", "DOMICILIO01_V0008"), df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
    d_sew <- extract_column_multi(df_household, NULL, "hh_sewage", c("V0016", "DOMICILIO01_V0016"), df_ref_hh = p_hh$total_hh, min_prop = 0, max_prop = 1)
  }
  
  df_std <- p_pop |>
    dplyr::left_join(p_hh,  by = "CODE_TRACT") |>
    dplyr::left_join(p_inc, by = "CODE_TRACT") |>
    dplyr::left_join(d_apt, by = "CODE_TRACT") |>
    dplyr::left_join(d_sew, by = "CODE_TRACT") |>
    dplyr::mutate(
      census_year = year,
      total_pop   = as.numeric(total_pop),
      total_hh    = as.numeric(total_hh),
      avg_income  = as.numeric(avg_income),
      hh_apt      = as.numeric(hh_apt),
      hh_sewage   = as.numeric(hh_sewage)
    )
  
  sf_year <- sf_tracts |>
    dplyr::left_join(df_std, by = c("code_tract" = "CODE_TRACT")) |>
    dplyr::filter(!is.na(total_pop) & total_pop > 0)
  
  if (nrow(sf_year) == 0) return(NULL)
  
  # Imputation based on geographic neighbors
  centroids <- suppressWarnings(sf::st_centroid(sf_year))
  knn_obj   <- spdep::knearneigh(sf::st_coordinates(centroids), k = 5)
  nb_obj    <- spdep::knn2nb(knn_obj)
  
  sf_year$avg_income <- spatial_impute(sf_year$avg_income, nb_obj)
  sf_year$hh_apt     <- spatial_impute(sf_year$hh_apt,     nb_obj)
  sf_year$hh_sewage  <- spatial_impute(sf_year$hh_sewage,  nb_obj)
  
  sf_year <- sf_year |>
    dplyr::mutate(
      log_income     = log(avg_income + 1),
      prop_apt       = ifelse(total_hh > 0, pmin((hh_apt    / total_hh) * 100, 100), NA_real_),
      prop_sewage    = ifelse(total_hh > 0, pmin((hh_sewage / total_hh) * 100, 100), NA_real_),
      persons_per_hh = ifelse(total_hh > 0, total_pop / total_hh, NA_real_)
    )
  
  # Outlier treatment (Winsorization at 99th percentile)
  q99_pph <- quantile(sf_year$persons_per_hh, 0.99, na.rm = TRUE)
  sf_year$persons_per_hh <- ifelse(sf_year$persons_per_hh > q99_pph, q99_pph, sf_year$persons_per_hh)
  
  sf_year <- sf_year |>
    dplyr::mutate(
      area_sqkm    = as.numeric(sf::st_area(sf_year)) / 1e6,
      density_sqkm = total_pop / area_sqkm,
      log_density  = log(density_sqkm + 1)
    )
  
  # Accessibility and informality variables
  sf_year$dist_train_km          <- as.numeric(sf::st_distance(centroids, sf_train[sf::st_nearest_feature(centroids, sf_train), ], by_element = TRUE)) / 1000
  sf_year$dist_bus_km            <- as.numeric(sf::st_distance(centroids, sf_bus[sf::st_nearest_feature(centroids, sf_bus), ],      by_element = TRUE)) / 1000
  sf_year$pct_area_slum          <- calc_informality_overlap(sf_year, sf_slum)
  sf_year$pct_area_tenement      <- calc_informality_overlap(sf_year, sf_tenement)
  sf_year$pct_area_irregular_lot <- calc_informality_overlap(sf_year, sf_irregular_lot)
  
  return(sf_year)
}

#-------------------------------------------------------------------------------
# ALGORITHM UTILITIES
#-------------------------------------------------------------------------------

calc_adaptive_bisquare <- function(d_s, d_t, h_s, h_t) {
  w_s <- ifelse(d_s <= h_s, (1 - (d_s / h_s)^2)^2, 0)
  w_t <- ifelse(d_t <= h_t, (1 - (d_t / h_t)^2)^2, 0)
  return(w_s * w_t)
}

golden_section <- function(f, lower, upper, tol = 1e-3, maxiter = 50) {
  gr <- (sqrt(5) - 1) / 2
  a <- lower; b <- upper
  c <- b - gr * (b - a); d <- a + gr * (b - a)
  fc <- f(c); fd <- f(d)
  iter <- 0
  while ((b - a) > tol && iter < maxiter) {
    if (fc < fd) {
      b <- d; d <- c; fd <- fc
      c <- b - gr * (b - a); fc <- f(c)
    } else {
      a <- c; c <- d; fc <- fd
      d <- a + gr * (b - a); fd <- f(d)
    }
    iter <- iter + 1
  }
  return((a + b) / 2)
}

validate_model_input <- function(data_sf, formula, time_var, required_cols = NULL) {
  if (!inherits(data_sf, "sf")) stop("data_sf must be an sf object.")
  df <- sf::st_drop_geometry(data_sf)
  if (!inherits(formula, "formula")) stop("formula must be an R formula.")
  
  y_var  <- as.character(formula[[2]])
  x_vars <- attr(terms(formula), "term.labels")
  needed  <- c(y_var, x_vars, time_var, required_cols)
  missing <- setdiff(needed, names(df))
  
  if (length(missing) > 0) stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  
  numeric_needed <- c(y_var, x_vars)
  non_numeric    <- numeric_needed[!sapply(df[numeric_needed], is.numeric)]
  
  if (length(non_numeric) > 0) stop(sprintf("Columns must be numeric: %s", paste(non_numeric, collapse = ", ")))
  TRUE
}

compute_local_t <- function(df, var_name) {
  beta_col <- paste0("Beta_", var_name)
  se_col   <- paste0("SE_",   var_name)
  if (se_col %in% names(df)) {
    return(df[[beta_col]] / df[[se_col]])
  } else {
    return(df[[beta_col]] / sd(df[[beta_col]], na.rm = TRUE))
  }
}

#-------------------------------------------------------------------------------
# BENCHMARK MGTWR ALGORITHMS 
#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 2SCALL
# ------------------------------------------------------------------------------
mgtwr_2scall_impl <- function(formula, data_sf, time_var, c_s = 0.8, c_t = 0.8) {
  x_vars <- attr(terms(formula), "term.labels")
  y_var  <- as.character(formula[[2]])
  df     <- sf::st_drop_geometry(data_sf)
  n      <- nrow(df)
  coords <- sf::st_coordinates(sf::st_centroid(data_sf))
  times  <- as.numeric(as.character(df[[time_var]]))
  Y      <- as.numeric(df[[y_var]])
  X      <- cbind(1, as.matrix(df[, x_vars, drop = FALSE]))
  P      <- ncol(X)
  
  nn_global <- RANN::nn2(coords, coords, k = min(n, 150))
  t_range <- max(times) - min(times); if (t_range == 0) t_range <- 1
  
  h_t_g <- t_range * 0.8 * c_t 
  if (h_t_g == 0) h_t_g <- 1
  
  beta_step1 <- matrix(0, nrow = n, ncol = P)
  for (i in seq_len(n)) {
    idx <- nn_global$nn.idx[i, ]
    h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * c_s, 1e-4)
    
    w_v <- calc_adaptive_bisquare(nn_global$nn.dists[i, ], abs(times[idx] - times[i]), h_s_i, h_t_g)
    v <- which(w_v > 1e-4)
    
    if (length(v) < 4 * P) {
      v <- 1:min(length(idx), 4 * P)
      w_v[v[w_v[v] <= 1e-4]] <- 1e-4
    }
    idx_v <- idx[v]; w_v <- w_v[v]
    
    du <- (coords[idx_v, 1] - coords[i, 1]) / h_s_i
    dv <- (coords[idx_v, 2] - coords[i, 2]) / h_s_i
    dt <- (times[idx_v] - times[i]) / h_t_g
    X_sub <- X[idx_v, , drop = FALSE]
    X_0 <- cbind(X_sub, X_sub * du, X_sub * dv, X_sub * dt)
    
    M <- crossprod(X_0, X_0 * w_v) + diag(1e-5, ncol(X_0))
    inv <- tryCatch(chol2inv(chol(M)), error = function(e) MASS::ginv(M))
    beta_step1[i, ] <- (inv %*% crossprod(X_0, w_v * Y[idx_v]))[1:P]
  }
  
  beta_step2 <- matrix(0, nrow = n, ncol = P)
  for (k in seq_len(P)) {
    Y_part <- Y - rowSums(beta_step1[, -k, drop = FALSE] * X[, -k, drop = FALSE])
    for (i in seq_len(n)) {
      idx <- nn_global$nn.idx[i, ]
      h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * c_s, 1e-4)
      
      w_v <- calc_adaptive_bisquare(nn_global$nn.dists[i, ], abs(times[idx] - times[i]), h_s_i / c_s, h_t_g / c_t)
      v <- which(w_v > 1e-4)
      if (length(v) < 4) {
        v <- 1:min(length(idx), 4)
        w_v[v[w_v[v] <= 1e-4]] <- 1e-4
      }
      idx_v <- idx[v]; w_v <- w_v[v]
      
      du <- (coords[idx_v, 1] - coords[i, 1]) / (h_s_i / c_s)
      dv <- (coords[idx_v, 2] - coords[i, 2]) / (h_s_i / c_s)
      dt <- (times[idx_v] - times[i]) / (h_t_g / c_t)
      
      X_k_0 <- cbind(X[idx_v, k], X[idx_v, k] * du, X[idx_v, k] * dv, X[idx_v, k] * dt)
      M <- crossprod(X_k_0, X_k_0 * w_v) + diag(1e-5, 4)
      inv_k <- tryCatch(chol2inv(chol(M)), error = function(e) MASS::ginv(M))
      beta_step2[i, k] <- (inv_k %*% crossprod(X_k_0, w_v * Y_part[idx_v]))[1]
    }
  }
  return(list(coefs = beta_step2, yhat = rowSums(beta_step2 * X), residual = Y - rowSums(beta_step2 * X)))
}

# ------------------------------------------------------------------------------
# TDS-MGTWR
# ------------------------------------------------------------------------------
mgtwr_tds_impl <- function(formula, data_sf, time_var, c_s = 0.8, c_t = 0.8) {
  x_vars <- attr(terms(formula), "term.labels")
  y_var <- as.character(formula[[2]])
  df <- sf::st_drop_geometry(data_sf); n <- nrow(df)
  coords <- sf::st_coordinates(sf::st_centroid(data_sf)); times <- as.numeric(as.character(df[[time_var]]))
  Y <- as.numeric(df[[y_var]]); X <- cbind(1, as.matrix(df[, x_vars, drop = FALSE])); P <- ncol(X)
  
  nn_global <- RANN::nn2(coords, coords, k = min(n, round(n * 0.5)))
  t_range <- max(times) - min(times); if (t_range == 0) t_range <- 1
  
  h_t_init <- t_range * 0.8 * c_t 
  if (h_t_init == 0) h_t_init <- 1
  
  global_mod <- lm(formula, data = df)
  top_down_order <- order(abs(summary(global_mod)$coefficients[, "t value"]), decreasing = TRUE)
  
  beta_mat <- matrix(0, nrow = n, ncol = P); bws_mat <- matrix(NA, nrow = P, ncol = 2)
  
  for (i in seq_len(n)) {
    idx <- nn_global$nn.idx[i, ]
    h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * c_s, 1e-4) 
    
    w_v <- calc_adaptive_bisquare(nn_global$nn.dists[i, ], abs(times[idx] - times[i]), h_s_i, h_t_init)
    v <- which(w_v > 1e-5) 
    
    if (length(v) < P) {
      v <- 1:min(length(idx), P + 2)
      w_v[v[w_v[v] <= 1e-5]] <- 1e-5 
    }
    
    idx_v <- idx[v]; w_v <- w_v[v]; X_sub <- X[idx_v, , drop = FALSE]
    M <- crossprod(X_sub, X_sub * w_v) + diag(1e-6, P)
    inv <- tryCatch(chol2inv(chol(M)), error = function(e) MASS::ginv(M))
    beta_mat[i, ] <- (inv %*% crossprod(X_sub, w_v * Y[idx_v]))[1:P]
  }
  
  eval_tds_aicc <- function(hs_mult, ht_try, k_idx, Y_part) {
    err_sq_sum <- 0; trace_S <- 0
    for (i in seq_len(n)) {
      idx <- nn_global$nn.idx[i, ]
      h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * hs_mult, 1e-4)
      
      w_v <- calc_adaptive_bisquare(nn_global$nn.dists[i, ], abs(times[idx] - times[i]), h_s_i, ht_try)
      v <- which(w_v > 1e-5)
      
      if (length(v) < 3) {
        v <- 1:min(length(idx), 3)
        w_v[v[w_v[v] <= 1e-5]] <- 1e-5
      }
      
      idx_v <- idx[v]; w_v <- w_v[v]; x_k <- X[idx_v, k_idx]
      inv_k <- 1 / (sum(w_v * x_k^2) + 1e-6)
      y_hat_local <- (inv_k * sum(x_k * w_v * Y_part[idx_v])) * X[i, k_idx]
      err_sq_sum <- err_sq_sum + (Y_part[i] - y_hat_local)^2
      trace_S <- trace_S + (w_v[1] * (X[i, k_idx]^2) * inv_k)
    }
    sigma2 <- err_sq_sum / n
    if (sigma2 <= 0 || trace_S >= n - 1) return(Inf)
    return(n * log(sigma2) + n * log(2 * pi) + n + (2 * trace_S * (trace_S + 1)) / (n - trace_S - 1))
  }
  
  for (order_idx in seq_along(top_down_order)) {
    k <- top_down_order[order_idx]
    Y_part <- Y - rowSums(beta_mat[, -k, drop = FALSE] * X[, -k, drop = FALSE])
    
    hs_mult_opt <- golden_section(function(mult) eval_tds_aicc(mult, h_t_init, k, Y_part), lower = 0.1, upper = 10)
    ht_opt <- golden_section(function(ht) eval_tds_aicc(hs_mult_opt, ht, k, Y_part), lower = t_range * 0.1, upper = t_range * 2)
    bws_mat[k, ] <- c(hs_mult_opt, ht_opt)
    
    for (i in seq_len(n)) {
      idx <- nn_global$nn.idx[i, ]
      h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * hs_mult_opt, 1e-4)
      
      w_v <- calc_adaptive_bisquare(nn_global$nn.dists[i, ], abs(times[idx] - times[i]), h_s_i, ht_opt)
      v <- which(w_v > 1e-5)
      
      if (length(v) < 3) {
        v <- 1:min(length(idx), 3)
        w_v[v[w_v[v] <= 1e-5]] <- 1e-5
      }
      
      idx_v <- idx[v]; w_v <- w_v[v]; x_k <- X[idx_v, k]
      inv_k <- 1 / (sum(w_v * x_k^2) + 1e-6)
      beta_mat[i, k] <- inv_k * sum(x_k * w_v * Y_part[idx_v])
    }
  }
  yhat <- rowSums(beta_mat * X)
  return(list(coefs = beta_mat, yhat = yhat, residual = Y - yhat, optimal_bws = bws_mat))
}

# ------------------------------------------------------------------------------
# Back-fitting
# ------------------------------------------------------------------------------
mgtwr_fotheringham_impl <- function(formula, data_sf, time_var, c_s = 0.8, c_t = 0.8) {
  x_vars <- attr(terms(formula), "term.labels")
  y_var  <- as.character(formula[[2]])
  df     <- sf::st_drop_geometry(data_sf)
  n      <- nrow(df)
  coords <- sf::st_coordinates(sf::st_centroid(data_sf))
  times  <- as.numeric(as.character(df[[time_var]]))
  Y      <- as.numeric(df[[y_var]])
  X      <- cbind(1, as.matrix(df[, x_vars, drop = FALSE]))
  P      <- ncol(X)
  colnames(X)[1] <- "Intercept"
  
  nn_global <- RANN::nn2(coords, coords, k = min(n, round(n * 0.5)))
  t_range <- max(times) - min(times)
  if (t_range == 0) t_range <- 1
  
  h_t_init <- t_range * 0.8 * c_t
  if (h_t_init == 0) h_t_init <- 1
  
  beta_mat <- matrix(0, nrow = n, ncol = P)
  
  for (i in seq_len(n)) {
    idx <- nn_global$nn.idx[i, ]
    h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * c_s, 1e-4)
    d_s <- nn_global$nn.dists[i, ]
    d_t <- abs(times[idx] - times[i])
    
    w_v <- calc_adaptive_bisquare(d_s, d_t, h_s_i, h_t_init)
    v <- which(w_v > 1e-5) 
    
    if (length(v) < P) {
      v <- 1:min(length(idx), P + 2)
      w_v[v[w_v[v] <= 1e-5]] <- 1e-5 
    }
    
    idx_v <- idx[v]; w_v <- w_v[v]; X_sub <- X[idx_v, , drop = FALSE]
    M <- crossprod(X_sub, X_sub * w_v) + diag(1e-6, P)
    inv <- tryCatch(chol2inv(chol(M)), error = function(e) MASS::ginv(M))
    beta_mat[i, ] <- (inv %*% crossprod(X_sub, w_v * Y[idx_v]))[1:P]
  }
  
  bws_mat <- matrix(NA, nrow = P, ncol = 2)
  rownames(bws_mat) <- colnames(X); colnames(bws_mat) <- c("h_s_opt", "h_t_opt")
  
  eval_mgtwr_aicc <- function(hs_mult, ht_try, k_idx, Y_part) {
    err_sq_sum <- 0; trace_S <- 0
    for (i in seq_len(n)) {
      idx <- nn_global$nn.idx[i, ]
      h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * hs_mult, 1e-4)
      d_s <- nn_global$nn.dists[i, ]
      d_t <- abs(times[idx] - times[i])
      
      w_v <- calc_adaptive_bisquare(d_s, d_t, h_s_i, ht_try)
      v <- which(w_v > 1e-5)
      
      if (length(v) < 3) {
        v <- 1:min(length(idx), 3)
        w_v[v[w_v[v] <= 1e-5]] <- 1e-5
      }
      
      idx_v <- idx[v]; w_v <- w_v[v]; x_k <- X[idx_v, k_idx]
      inv_k <- 1 / (sum(w_v * x_k^2) + 1e-6)
      y_hat_local <- (inv_k * sum(x_k * w_v * Y_part[idx_v])) * X[i, k_idx]
      err_sq_sum <- err_sq_sum + (Y_part[i] - y_hat_local)^2
      trace_S <- trace_S + (w_v[1] * (X[i, k_idx]^2) * inv_k)
    }
    sigma2 <- err_sq_sum / n
    if (sigma2 <= 0 || trace_S >= n - 1) return(Inf)
    return(n * log(sigma2) + n * log(2 * pi) + n + (2 * trace_S * (trace_S + 1)) / (n - trace_S - 1))
  }
  
  max_iter <- 50
  tol <- 1e-5
  soc_f <- Inf
  iter <- 0
  
  se_mat <- matrix(0, nrow = n, ncol = P)
  
  while (soc_f > tol && iter < max_iter) {
    iter <- iter + 1
    beta_old <- beta_mat
    
    for (k in seq_len(P)) {
      Y_part <- Y - rowSums(beta_mat[, -k, drop = FALSE] * X[, -k, drop = FALSE])
      
      hs_mult_opt <- golden_section(function(mult) eval_mgtwr_aicc(mult, h_t_init, k, Y_part), lower = 0.1, upper = 10)
      ht_opt      <- golden_section(function(ht) eval_mgtwr_aicc(hs_mult_opt, ht, k, Y_part), lower = t_range * 0.1, upper = t_range * 2)
      
      bws_mat[k, ] <- c(hs_mult_opt, ht_opt)
      
      for (i in seq_len(n)) {
        idx <- nn_global$nn.idx[i, ]
        h_s_i <- max(nn_global$nn.dists[i, min(50, ncol(nn_global$nn.dists))] * hs_mult_opt, 1e-4)
        d_s <- nn_global$nn.dists[i, ]
        d_t <- abs(times[idx] - times[i])
        
        w_v <- calc_adaptive_bisquare(d_s, d_t, h_s_i, ht_opt)
        v <- which(w_v > 1e-5)
        
        if (length(v) < 3) {
          v <- 1:min(length(idx), 3)
          w_v[v[w_v[v] <= 1e-5]] <- 1e-5
        }
        
        idx_v <- idx[v]; w_v <- w_v[v]; x_k <- X[idx_v, k]
        
        M <- sum(w_v * x_k^2) + 1e-6
        inv_k <- 1 / M
        beta_mat[i, k] <- inv_k * sum(x_k * w_v * Y_part[idx_v])
        
        if (soc_f <= tol || iter == max_iter) {
          tr_S <- w_v[1] * (X[i, k]^2) * inv_k 
          local_var <- sum(w_v * (Y_part[idx_v] - (x_k * beta_mat[i, k]))^2) / max(1, sum(w_v) - tr_S)
          se_mat[i, k] <- sqrt(max(0, local_var * inv_k))
        }
      }
    }
    
    num <- sum((beta_mat * X - beta_old * X)^2)
    den <- sum((beta_mat * X)^2)
    soc_f <- sqrt(num / (den + 1e-12))
  }
  
  yhat <- rowSums(beta_mat * X)
  return(list(coefs = beta_mat, se = se_mat, yhat = yhat, residual = Y - yhat, optimal_bws = bws_mat, iterations = iter))
}

#-------------------------------------------------------------------------------
# DIAGNOSTICS & PLOTTING UTILITIES
#-------------------------------------------------------------------------------

watson_U2 <- function(u) {
  u <- sort(u[!is.na(u)]); n <- length(u)
  W2 <- (1/(12*n)) + sum((u - (2*(1:n)-1)/(2*n))^2)
  return(list(U2 = W2 - n*(mean(u) - 0.5)^2, W2 = W2))
}

neyman_ledwina <- function(u, Smax = 6) {
  u <- u[!is.na(u)]; n <- length(u)
  Xbase <- function(k, x) sapply(0:k, function(p) (2*x - 1)^p)
  bic_vals <- numeric(Smax); stats <- numeric(Smax)
  for(S in 1:Smax) {
    X <- Xbase(S, u); Q <- qr.Q(qr(X)); B <- Q[, -1, drop = FALSE]
    Tj <- colSums(B) / sqrt(n); stats[S] <- sum(Tj^2)
    bic_vals[S] <- -stats[S] + S * log(n)
  }
  return(list(Sstar = which.min(bic_vals), max_stat = max(stats)))
}

tb2_entropy_proxy <- function(u, m = 3) {
  u_s <- sort(u[!is.na(u)]); n <- length(u_s)
  m <- min(m, floor((n-1)/2))
  diffs <- sapply((m+1):(n-m), function(i) u_s[i+m] - u_s[i-m])
  diffs[diffs <= 0] <- 1e-12
  return(list(TB2 = sum(log(diffs))))
}

compute_continuous_pit <- function(residuals) {
  sigma <- stats::mad(residuals, constant = 1.4826, na.rm = TRUE)
  pit <- pnorm(residuals, mean = 0, sd = sigma)
  pmax(pmin(pit, 0.999999), 0.000001) 
}

plot_pit_histogram <- function(pit_vals, model_name, u2_val) {
  df <- data.frame(PIT = pit_vals)
  ggplot2::ggplot(df, ggplot2::aes(x = PIT)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), 
                            bins = 30, fill = "slategray", color = "white", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", linewidth = 1) +
    ggplot2::labs(title = model_name, 
                  subtitle = sprintf("Watson U²: %.4f", u2_val), 
                  x = "Probability Integral Transform (PIT)", y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))
}

calc_model_metrics <- function(y, yhat, coords = NULL, k = 5, exec_time = NA) {
  n <- length(y); rss <- sum((y - yhat)^2, na.rm = TRUE); tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  moran_i <- NA_real_; moran_p <- NA_real_
  if (!is.null(coords)) {
    listw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k)), style = "W", zero.policy = TRUE)
    moran <- tryCatch(spdep::moran.test(y - yhat, listw, zero.policy = TRUE), error = function(e) NULL)
    if (!is.null(moran)) { moran_i <- moran$estimate[1]; moran_p <- moran$p.value }
  }
  return(tibble::tibble(RSS = rss, R2 = 1 - (rss / tss), RMSE = sqrt(rss / n), MAE = mean(abs(y - yhat), na.rm = TRUE), MoranI = moran_i, MoranP = moran_p, Time_sec = round(exec_time, 2)))
}

compute_worm_pit <- function(residuals) {
  n <- length(residuals)
  sigma <- stats::mad(residuals, constant = 1.4826, na.rm = TRUE)
  tibble::tibble(residual = residuals) |>
    dplyr::arrange(residual) |>
    dplyr::mutate(
      PIT           = pnorm(residual, mean = 0, sd = sigma),
      theoretical_p = (dplyr::row_number() - 0.5) / n,
      z_expected    = qnorm(theoretical_p),
      z_obs         = residual / sigma,
      worm_deviation = z_obs - z_expected
    )
}

create_wormplot <- function(pit_values, model_name) {
  qres <- qnorm(pit_values)
  qres <- sort(qres[!is.na(qres) & !is.infinite(qres)])
  n <- length(qres)
  p <- (1:n - 0.5) / n
  z <- qnorm(p)
  dz <- qres - z 
  
  se <- sqrt(p * (1 - p) / n)
  upper <- qnorm(pmin(p + 1.96 * se, 0.9999)) - z
  lower <- qnorm(pmax(p - 1.96 * se, 0.0001)) - z
  
  wp_df <- data.frame(z = z, dz = dz, upper = upper, lower = lower)
  
  ggplot2::ggplot(wp_df, ggplot2::aes(x = z, y = dz)) +
    ggplot2::geom_point(alpha = 0.6, color = "#2c7fb8") +
    ggplot2::geom_line(ggplot2::aes(y = upper), linetype = "dashed", color = "#e34a33", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = lower), linetype = "dashed", color = "#e34a33", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Worm Plot -", model_name), subtitle = "Model Calibration & Fit",
                  x = "Theoretical Normal Quantiles", y = "Deviation")
}

plot_sf_variable <- function(sf_obj, var, facet = NULL, districts = NULL, palette = "viridis", title = NULL, na.value = NA) {
  p <- ggplot2::ggplot(sf_obj) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA) +
    ggplot2::scale_fill_viridis_c(option = palette, na.value = na.value) +
    ggplot2::labs(title = title, fill = var)
  if (!is.null(districts)) p <- p + ggplot2::geom_sf(data = districts, fill = NA, color = "black", linewidth = 0.3)
  if (!is.null(facet))     p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet)))
  p + ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))
}

plot_masked_coefficients <- function(sf_obj, beta_var, districts = NULL, title = "Extreme Deviations from Global Mean") {
  t_stat <- sf_obj[[beta_var]] / sd(sf_obj[[beta_var]], na.rm = TRUE)
  sf_obj$sig_beta <- ifelse(abs(t_stat) >= 1.96, sf_obj[[beta_var]], NA)
  p <- ggplot2::ggplot(sf_obj) +
    ggplot2::geom_sf(fill = "gray90", color = NA) +
    ggplot2::geom_sf(ggplot2::aes(fill = sig_beta), color = NA) +
    ggplot2::scale_fill_viridis_c(option = "mako", na.value = "transparent") +
    ggplot2::facet_wrap(~census_year) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))
  if (!is.null(districts)) p <- p + ggplot2::geom_sf(data = districts, fill = NA, color = "black", linewidth = 0.2)
  return(p)
}

plot_unmasked <- function(model_data_real_out, districts_clean, var, title, option, name, file_name, out_dir = "outputs") {
  p <- ggplot2::ggplot(model_data_real_out) +
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 11, progress = "none") +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA, alpha = 0.8) +
    ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(option = option, direction = -1, name = name) +
    ggplot2::facet_wrap(~ census_year, ncol = 3) +
    ggplot2::labs(title = title, caption = unmasked_caption) + 
    academic_map_theme
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(out_dir, file_name), plot = p, width = 12, height = 5, dpi = 300, bg = "white")
}

plot_unmasked_custom <- function(model_data_real_out, districts_clean, var, title, option, name, file_name, out_dir) {
  p <- ggplot2::ggplot(model_data_real_out) +
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 11, progress = "none") +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA, alpha = 0.8) +
    ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(option = option, direction = -1, name = name) +
    ggplot2::facet_wrap(~ census_year, ncol = 3) +
    ggplot2::labs(title = title, caption = "Note: Standardized variables. Parameters estimated via TDS-MGTWR.") + 
    theme_publication() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  ggplot2::ggsave(file.path(out_dir, file_name), plot = p, width = 12, height = 5, dpi = 300, bg = "white")
}

plot_masked <- function(model_data_real_out, districts_clean, var, title, option, name, file_name, out_dir = "outputs") {
  p <- ggplot2::ggplot(model_data_real_out) +
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 11, progress = "none") +
    ggplot2::geom_sf(data = dplyr::filter(model_data_real_out, is.na(.data[[var]])), fill = "gray90", color = NA, alpha = 0.8) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA, alpha = 0.8) +
    ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(option = option, direction = -1, name = name, na.value = "transparent") +
    ggplot2::facet_wrap(~ census_year, ncol = 3) +
    ggplot2::labs(title = title, caption = masked_caption) + 
    academic_map_theme
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(out_dir, file_name), plot = p, width = 12, height = 5, dpi = 300, bg = "white")
}

plot_masked_custom <- function(model_data_real_out, districts_clean, var, title, option, name, file_name, out_dir) {
  p <- ggplot2::ggplot(model_data_real_out) +
    ggspatial::annotation_map_tile(type = "cartolight", zoom = 11, progress = "none") +
    ggplot2::geom_sf(data = dplyr::filter(model_data_real_out, is.na(.data[[var]])), fill = "gray90", color = NA, alpha = 0.8) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), color = NA, alpha = 0.8) +
    ggplot2::geom_sf(data = districts_clean, fill = NA, color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(option = option, direction = -1, name = name, na.value = "transparent") +
    ggplot2::facet_wrap(~ census_year, ncol = 3) +
    ggplot2::labs(title = title, caption = "Note: Gray areas represent coefficients lacking extreme spatial deviation (>1.96 SD).") + 
    theme_publication() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  ggplot2::ggsave(file.path(out_dir, file_name), plot = p, width = 12, height = 5, dpi = 300, bg = "white")
}

#-------------------------------------------------------------------------------
# MODEL WRAPPERS
#-------------------------------------------------------------------------------

run_mgtwr_2scall <- function(formula, data_sf, time_var, cache_path = NULL) {
  if (!is.null(cache_path) && file.exists(cache_path)) return(readRDS(cache_path))
  t0 <- Sys.time()
  res <- mgtwr_2scall_impl(formula, data_sf, time_var)
  res$exec_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(cache_path)) saveRDS(res, cache_path)
  return(res)
}

run_mgtwr_tds <- function(formula, data_sf, time_var, cache_path = NULL) {
  if (!is.null(cache_path) && file.exists(cache_path)) return(readRDS(cache_path))
  t0 <- Sys.time()
  res <- mgtwr_tds_impl(formula, data_sf, time_var)
  res$exec_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(cache_path)) saveRDS(res, cache_path)
  return(res)
}

run_mgtwr_fotheringham <- function(formula, data_sf, time_var, cache_path = NULL) {
  if (!is.null(cache_path) && file.exists(cache_path)) return(readRDS(cache_path))
  t0 <- Sys.time()
  res <- mgtwr_fotheringham_impl(formula, data_sf, time_var)
  res$exec_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(cache_path)) saveRDS(res, cache_path)
  return(res)
}

#-------------------------------------------------------------------------------
# EVALUATION
#-------------------------------------------------------------------------------

evaluate_and_plot_mgtwr <- function(model_data, districts_clean = NULL, prefix_name = "model", formula, time_var = "census_year", map_var = NULL, out_dir = "outputs", models_to_run = c("2SCALL", "TDS", "F2017")) {
  validate_model_input(model_data, formula, time_var)
  
  y_var <- as.character(formula[[2]])
  x_vars <- attr(terms(formula), "term.labels")
  coef_names <- c("Intercept", x_vars)
  coords <- sf::st_coordinates(sf::st_centroid(model_data))
  
  mod_2scall <- NULL; mod_tds <- NULL; mod_f2017 <- NULL
  metrics_list <- list()
  
  models_to_run <- toupper(models_to_run)
  
  if ("2SCALL" %in% models_to_run) {
    mod_2scall <- run_mgtwr_2scall(formula, model_data, time_var)
    df_coefs_2scall <- as.data.frame(mod_2scall$coefs); colnames(df_coefs_2scall) <- paste0("Beta_", coef_names, "_2scall")
    model_data <- dplyr::bind_cols(model_data, df_coefs_2scall)
    metrics_list[[length(metrics_list) + 1]] <- dplyr::mutate(calc_model_metrics(model_data[[y_var]], mod_2scall$yhat, coords, exec_time = mod_2scall$exec_time), Model = "MGTWR_2SCALL")
  }
  
  if ("TDS" %in% models_to_run) {
    mod_tds <- run_mgtwr_tds(formula, model_data, time_var)
    df_coefs_tds <- as.data.frame(mod_tds$coefs); colnames(df_coefs_tds) <- paste0("Beta_", coef_names, "_tds")
    model_data <- dplyr::bind_cols(model_data, df_coefs_tds)
    metrics_list[[length(metrics_list) + 1]] <- dplyr::mutate(calc_model_metrics(model_data[[y_var]], mod_tds$yhat, coords, exec_time = mod_tds$exec_time), Model = "MGTWR_TDS")
  }
  
  #F2017=BF
  if ("F2017" %in% models_to_run) {
    mod_f2017 <- run_mgtwr_fotheringham(formula, model_data, time_var)
    df_coefs_f2017 <- as.data.frame(mod_f2017$coefs); colnames(df_coefs_f2017) <- paste0("Beta_", coef_names)
    df_se_f2017    <- as.data.frame(mod_f2017$se);    colnames(df_se_f2017)    <- paste0("SE_", coef_names)
    model_data  <- dplyr::bind_cols(model_data, df_coefs_f2017, df_se_f2017)
    model_data$yhat_f2017 <- mod_f2017$yhat; model_data$resid_f2017 <- mod_f2017$residual
    metrics_list[[length(metrics_list) + 1]] <- dplyr::mutate(calc_model_metrics(model_data[[y_var]], mod_f2017$yhat, coords, exec_time = mod_f2017$exec_time), Model = "MGTWR_F2017")
  }
  
  metrics_df <- dplyr::bind_rows(metrics_list)
  if (nrow(metrics_df) > 0) {
    print(knitr::kable(metrics_df, caption = paste("Performance Comparison -", prefix_name)))
  }
  
  if ("F2017" %in% models_to_run) {
    sd_y <- stats::sd(sf::st_drop_geometry(model_data)[[y_var]], na.rm = TRUE)
    model_data[["Beta_Intercept_Original"]] <- model_data[["Beta_Intercept"]]
    
    for (var in x_vars) {
      beta_col <- paste0("Beta_", var)
      orig_col <- paste0(beta_col, "_Original")
      if (beta_col %in% names(model_data) && var %in% names(model_data)) {
        sd_x <- stats::sd(sf::st_drop_geometry(model_data)[[var]], na.rm = TRUE)
        model_data[[orig_col]] <- if (is.na(sd_x) || sd_x == 0) model_data[[beta_col]] else model_data[[beta_col]] * (sd_y / sd_x)
      }
    }
    
    coef_match <- grep("^Beta_.*_Original$", names(model_data), value = TRUE)
    if (length(coef_match) > 0) {
      coef_summary <- model_data |> sf::st_drop_geometry() |> dplyr::select(dplyr::all_of(coef_match)) |>
        tidyr::pivot_longer(dplyr::everything(), names_to = "Variable", values_to = "Coefficient") |>
        dplyr::group_by(Variable) |>
        dplyr::summarise(
          Min = min(Coefficient, na.rm = TRUE), `1stQu` = stats::quantile(Coefficient, 0.25, na.rm = TRUE),
          Median = stats::median(Coefficient, na.rm = TRUE), `3rdQu` = stats::quantile(Coefficient, 0.75, na.rm = TRUE),
          Max = max(Coefficient, na.rm = TRUE), SD = stats::sd(Coefficient, na.rm = TRUE)
        ) |> dplyr::arrange(Variable)
      
      if (requireNamespace("gt", quietly = TRUE)) {
        gt_tbl <- gt::gt(coef_summary) |> gt::fmt_number(columns = 2:7, decimals = 4)
        print(gt_tbl)
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
        gt::gtsave(gt_tbl, filename = file.path(out_dir, paste0("Table_Coefficients_", prefix_name, ".docx")))
      }
    }
    
    df_diag <- compute_worm_pit(model_data$resid_f2017)
    p_pit <- ggplot2::ggplot(df_diag, ggplot2::aes(x = PIT)) +
      ggplot2::geom_histogram(bins = 30, fill = "slategray") + ggplot2::labs(title = "PIT Histogram (F2017)")
    ggplot2::ggsave(filename = file.path(out_dir, paste0("Fig_PIT_", prefix_name, ".png")), plot = p_pit, width = 8, height = 4, bg = "white")
    
    qres <- stats::quantile(model_data$resid_f2017, probs = c(0.025, 0.975), na.rm = TRUE)
    model_data$Residual_Class <- factor(dplyr::case_when(
      model_data$resid_f2017 < qres[1] ~ "Over-prediction",
      model_data$resid_f2017 > qres[2] ~ "Under-prediction",
      TRUE ~ "Normal Fit"
    ), levels = c("Over-prediction", "Normal Fit", "Under-prediction"))
    
    p_out <- ggplot2::ggplot()
    if (!is.null(districts_clean)) p_out <- p_out + ggplot2::geom_sf(data = districts_clean, fill = NA, color = "gray50")
    p_out <- p_out + ggplot2::geom_sf(data = subset(model_data, Residual_Class != "Normal Fit"), ggplot2::aes(fill = Residual_Class), alpha = 0.9) +
      ggplot2::facet_wrap(as.formula(paste("~", time_var))) + ggplot2::labs(title = "Extreme Residuals (F2017)") + ggplot2::theme_void() + ggplot2::theme(legend.position = "bottom")
    ggplot2::ggsave(filename = file.path(out_dir, paste0("Fig_ExtremeResiduals_", prefix_name, ".png")), plot = p_out, width = 12, height = 6, bg = "white")
  }
  
  target_map_vars <- if (is.null(map_var)) x_vars else map_var
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (t_var in target_map_vars) {
    fill_col_2scall <- paste0("Beta_", t_var, "_2scall")
    fill_col_tds    <- paste0("Beta_", t_var, "_tds")
    fill_col_f2017  <- paste0("Beta_", t_var)
    
    p1 <- NULL; p2 <- NULL; p3 <- NULL
    if (fill_col_2scall %in% names(model_data)) p1 <- plot_sf_variable(model_data, fill_col_2scall, facet = time_var, districts = districts_clean, title = paste("Local Effect:", t_var, "(2SCALL)"))
    if (fill_col_tds %in% names(model_data))    p2 <- plot_sf_variable(model_data, fill_col_tds, facet = time_var, districts = districts_clean, title = paste("Local Effect:", t_var, "(TDS-MGTWR)"))
    if (fill_col_f2017 %in% names(model_data))  p3 <- plot_sf_variable(model_data, fill_col_f2017, facet = time_var, districts = districts_clean, title = paste("Local Effect:", t_var, "(Classic F2017)"))
    
    if (!is.null(p1) && !is.null(p2) && !is.null(p3)) {
      comp <- p1 / p2 / p3 + patchwork::plot_annotation(title = paste("Methodological Comparison -", prefix_name, "-", t_var))
      ggplot2::ggsave(filename = file.path(out_dir, paste0("Fig_Comparison_", prefix_name, "_", t_var, ".png")), plot = comp, width = 14, height = 15, bg = "white")
    }
  }
  
  return(list(mod_2scall = mod_2scall, mod_tds = mod_tds, mod_f2017 = mod_f2017, model_data = model_data, metrics = metrics_df))
}

#-------------------------------------------------------------------------------
# MONTE CARLO DATA GENERATING PROCESS (DGP)
#-------------------------------------------------------------------------------

generate_dgp <- function(seed, grid_size) {
  set.seed(seed)
  nx <- grid_size; ny <- grid_size 
  grid_sim <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(c(xmin=0, ymin=0, xmax=grid_size, ymax=grid_size))), n=c(nx, ny))
  sf_sim_base <- sf::st_sf(tract_id = 1:length(grid_sim), geometry = grid_sim)
  coords <- safe_centroids(sf_sim_base, method = "point_on_surface")
  
  panel_years <- c(1, 2, 3) 
  
  sim_list <- lapply(panel_years, function(t) {
    df <- sf_sim_base
    df$time <- t
    u <- coords[, 1]; v <- coords[, 2]
    
    df$X1 <- rnorm(nrow(df), mean = 0, sd = 1)
    df$X2 <- rnorm(nrow(df), mean = 0, sd = 1)
    
    df$True_Beta_Intercept <- 2 + 0.2*u + 0.2*v + 0.5*t 
    df$True_Beta_X1 <- sin(u / 1.5) + cos(v / 1.5) + 0.2*t
    df$True_Beta_X2 <- rep(1.5, nrow(df)) 
    
    epsilon <- rnorm(nrow(df), mean = 0, sd = 0.5) 
    df$Y <- df$True_Beta_Intercept + (df$True_Beta_X1 * df$X1) + (df$True_Beta_X2 * df$X2) + epsilon
    return(df)
  })
  do.call(rbind, sim_list) |> sf::st_as_sf()
}

#-------------------------------------------------------------------------------
# AUTOCORRELATION ANALYSIS
#-------------------------------------------------------------------------------
calculate_lisa_by_year <- function(yr, sf_data) {
  sf_yr <- sf_data |> filter(census_year == yr)
  nb <- poly2nb(sf_yr, queen = TRUE, snap = 0.001)
  
  no_neighbors <- which(card(nb) == 0)
  if(length(no_neighbors) > 0) {
    coords <- suppressWarnings(st_centroid(st_geometry(sf_yr)))
    for (i in no_neighbors) {
      distances <- st_distance(coords[i], coords)
      nb[[i]] <- order(as.numeric(distances))[2] 
    }
  }
  
  listw <- nb2listw(nb, style = 'W', zero.policy = TRUE)
  g_moran <- moran.test(sf_yr$log_density, listw = listw, randomisation = TRUE, zero.policy = TRUE)
  
  lisa <- localmoran(sf_yr$log_density, listw = listw, zero.policy = TRUE, alternative = "two.sided")
  z_var <- scale(sf_yr$log_density)[,1]
  z_lag <- lag.listw(listw, z_var, zero.policy = TRUE)
  
  sf_yr <- sf_yr |> mutate(
    z_var_col = z_var,
    z_lag_col = z_lag,
    pvalue = lisa[, 5], 
    cluster = case_when(
      z_var_col >  0 & z_lag_col >  0 & pvalue < 0.05 ~ "High-High (HH)",
      z_var_col <  0 & z_lag_col <  0 & pvalue < 0.05 ~ "Low-Low (LL)",
      z_var_col >  0 & z_lag_col <  0 & pvalue < 0.05 ~ "High-Low (HL)",
      z_var_col <  0 & z_lag_col >  0 & pvalue < 0.05 ~ "Low-High (LH)",
      TRUE                                            ~ "Not Significant"
    ),
    moran_global_i = g_moran$estimate[1],
    moran_p_value = g_moran$p.value
  ) |> 
    mutate(cluster = factor(cluster, levels = c("High-High (HH)", "Low-Low (LL)", "High-Low (HL)", "Low-High (LH)", "Not Significant")))
  
  return(sf_yr)
}

theme_n <- function() {
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey30"),
      axis.title       = ggplot2::element_text(face = "bold", size = 12),
      axis.text        = ggplot2::element_text(size = 10, color = "black"),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold"),
      legend.text      = ggplot2::element_text(size = 11),
      strip.text       = ggplot2::element_text(face = "bold", size = 12, color = "white"),
      strip.background = ggplot2::element_rect(fill = "#34495e", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
}