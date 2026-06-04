# Gaussian fitting functions for MAPLE
#
# These functions are adapted from the Gaussian fitting workflow implemented
# in PrInCE and extended for MAPLE analysis.
#
# MAPLE-specific extensions include:
#   - R2-based fit filtering
#   - RMSE-based fit filtering
#   - flexible pass rules using either or both R2 and RMSE
#   - optimisation by RMSE or R2
#   - random or guessed initialisation
#   - Gaussian component filtering by height, centre and variance
#   - early stopping once a fit passes the user-defined quality criteria
#
# Required packages:
#   PrInCE
#   purrr
#   progress

library(PrInCE)
library(purrr)
library(progress)


fit_gaussians_maple <- function(
    chromatogram,
    n_gaussians,
    max_iterations = 10,
    min_R_squared = 0.5,
    use_R_squared = TRUE,
    use_RMSE = FALSE,
    max_RMSE = Inf,
    pass_rule = c("both", "either"),
    optimize_by = c("RMSE", "R2"),
    method = c("guess", "random"),
    filter_gaussians_center = TRUE,
    filter_gaussians_height = 0.15,
    filter_gaussians_variance_min = 0.1,
    filter_gaussians_variance_max = 50
) {
  
  pass_rule <- match.arg(pass_rule)
  optimize_by <- match.arg(optimize_by)
  method <- match.arg(method)
  
  indices <- get_fraction_index(chromatogram)
  iter <- 0
  
  bestR2 <- -Inf
  bestRMSE <- Inf
  bestCoefs <- NULL
  bestPass <- FALSE
  
  while (iter < max_iterations && !bestPass) {
    iter <- iter + 1
    
    initial_conditions <- make_initial_conditions(
      chromatogram,
      n_gaussians,
      method
    )
    
    A <- initial_conditions$A
    mu <- initial_conditions$mu
    sigma <- initial_conditions$sigma
    
    p_model <- function(x, A, mu, sigma) {
      rowSums(sapply(seq_len(n_gaussians), function(i) {
        A[i] * exp(-((x - mu[i]) / sigma[i])^2)
      }))
    }
    
    fit <- tryCatch({
      suppressWarnings(
        nls(
          chromatogram ~ p_model(indices, A, mu, sigma),
          start = list(A = A, mu = mu, sigma = sigma),
          trace = FALSE,
          control = list(
            warnOnly = TRUE,
            minFactor = 1 / 2048
          )
        )
      )
    }, error = function(e) {
      e
    }, simpleError = function(e) {
      e
    })
    
    if ("error" %in% class(fit)) next
    
    coefs <- coef(fit)
    coefs <- split(coefs, rep(seq_len(3), each = n_gaussians))
    coefs <- setNames(coefs, c("A", "mu", "sigma"))
    
    if (filter_gaussians_variance_min > 0) {
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas < filter_gaussians_variance_min)
      if (length(drop) > 0) {
        coefs <- lapply(coefs, `[`, -drop)
      }
    }
    
    if (filter_gaussians_variance_max > 0) {
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas > filter_gaussians_variance_max)
      if (length(drop) > 0) {
        coefs <- lapply(coefs, `[`, -drop)
      }
    }
    
    if (filter_gaussians_center) {
      means <- coefs[["mu"]]
      drop <- which(means < 0 | means > length(chromatogram))
      if (length(drop) > 0) {
        coefs <- lapply(coefs, `[`, -drop)
      }
    }
    
    if (filter_gaussians_height > 0) {
      minHeight <- max(chromatogram) * filter_gaussians_height
      heights <- coefs[["A"]]
      drop <- which(heights < minHeight)
      if (length(drop) > 0) {
        coefs <- lapply(coefs, `[`, -drop)
      }
    }
    
    if (length(first(coefs)) == 0) next
    
    curveFit <- fit_curve(coefs, indices)
    
    R2 <- cor(chromatogram, curveFit)^2
    RMSE <- sqrt(mean((chromatogram - curveFit)^2, na.rm = TRUE))
    
    pass_R2 <- (!use_R_squared) || (R2 >= min_R_squared)
    pass_RMSE <- (!use_RMSE) || (RMSE <= max_RMSE)
    
    if (pass_rule == "both") {
      pass_current <- pass_R2 && pass_RMSE
    } else {
      pass_current <- pass_R2 || pass_RMSE
    }
    
    if (!pass_current) next
    
    if (optimize_by == "RMSE") {
      if (RMSE < bestRMSE) {
        bestRMSE <- RMSE
        bestR2 <- R2
        bestCoefs <- coefs
        bestPass <- TRUE
      }
    } else {
      if (R2 > bestR2) {
        bestR2 <- R2
        bestRMSE <- RMSE
        bestCoefs <- coefs
        bestPass <- TRUE
      }
    }
  }
  
  if (!is.null(bestCoefs)) {
    curveFit <- fit_curve(bestCoefs, indices)
  } else {
    curveFit <- NULL
    bestR2 <- NA_real_
    bestRMSE <- NA_real_
  }
  
  cols <- names(chromatogram)
  frac <- as.numeric(sub("^F", "", cols))
  offset <- min(frac) - 1
  
  if (!is.null(bestCoefs)) {
    bestCoefs$mu <- bestCoefs$mu + offset
  }
  
  results <- list(
    n_gaussians = n_gaussians,
    R2 = bestR2,
    RMSE = bestRMSE,
    iterations = iter,
    coefs = bestCoefs,
    curveFit = curveFit
  )
  
  return(results)
}


choose_gaussians_maple <- function(
    chromatogram,
    points = NULL,
    max_gaussians = 5,
    criterion = c("AICc", "AIC", "BIC"),
    max_iterations = 10,
    min_R_squared = 0.5,
    method = c("guess", "random"),
    filter_gaussians_center = TRUE,
    use_R_squared = TRUE,
    use_RMSE = FALSE,
    max_RMSE = Inf,
    pass_rule = c("both", "either"),
    optimize_by = c("RMSE", "R2"),
    filter_gaussians_height = 0.15,
    filter_gaussians_variance_min = 0.1,
    filter_gaussians_variance_max = 50
) {
  
  criterion <- match.arg(criterion)
  
  if (!is.null(points)) {
    max_gaussians <- min(max_gaussians, floor(points / 3))
  }
  
  fits <- list()
  
  for (n_gaussians in seq_len(max_gaussians)) {
    fits[[n_gaussians]] <- fit_gaussians_maple(
      chromatogram = chromatogram,
      n_gaussians = n_gaussians,
      max_iterations = max_iterations,
      min_R_squared = min_R_squared,
      use_R_squared = use_R_squared,
      use_RMSE = use_RMSE,
      max_RMSE = max_RMSE,
      pass_rule = pass_rule,
      optimize_by = optimize_by,
      method = method,
      filter_gaussians_center = filter_gaussians_center,
      filter_gaussians_height = filter_gaussians_height,
      filter_gaussians_variance_min = filter_gaussians_variance_min,
      filter_gaussians_variance_max = filter_gaussians_variance_max
    )
  }
  
  models <- map(fits, "coefs")
  drop <- map_lgl(models, is.null)
  fits <- fits[!drop]
  
  if (length(fits) == 0) {
    return(NULL)
  }
  
  coefs <- map(fits, "coefs")
  
  if (criterion == "AICc") {
    criteria <- lapply(coefs, gaussian_aicc, chromatogram)
  } else if (criterion == "AIC") {
    criteria <- lapply(coefs, gaussian_aic, chromatogram)
  } else if (criterion == "BIC") {
    criteria <- lapply(coefs, gaussian_bic, chromatogram)
  }
  
  best <- which.min(criteria)
  
  if (length(best) == 0) {
    return(NULL)
  } else {
    return(fits[[best]])
  }
}


build_gaussians_maple <- function(
    profile_matrix,
    min_points = 1,
    min_consecutive = 5,
    impute_NA = TRUE,
    smooth = TRUE,
    smooth_width = 4,
    max_gaussians = 5,
    criterion = c("AICc", "AIC", "BIC"),
    max_iterations = 50,
    min_R_squared = 0.5,
    use_R_squared = TRUE,
    use_RMSE = FALSE,
    max_RMSE = Inf,
    pass_rule = c("both", "either"),
    optimize_by = c("RMSE", "R2"),
    method = c("guess", "random"),
    filter_gaussians_center = TRUE,
    filter_gaussians_height = 0.15,
    filter_gaussians_variance_min = 0.5,
    filter_gaussians_variance_max = 50,
    show_progress = TRUE
) {
  
  if (is(profile_matrix, "MSnSet")) {
    profile_matrix <- exprs(profile_matrix)
  }
  
  filtered <- filter_profiles(
    profile_matrix,
    min_points = min_points,
    min_consecutive = min_consecutive
  )
  
  cleaned <- clean_profiles(
    filtered,
    impute_NA = impute_NA,
    smooth = smooth,
    smooth_width = smooth_width
  )
  
  cleaned <- cleaned[
    !(rownames(cleaned) == "" | is.na(rownames(cleaned))),
    , drop = FALSE
  ]
  
  gaussians <- list()
  proteins <- rownames(cleaned)
  P <- length(proteins)
  
  if (show_progress) {
    message(".. fitting Gaussian mixture models to ", P, " profiles")
    pb <- progress_bar$new(
      format = "fitting :what [:bar] :percent eta: :eta",
      clear = FALSE,
      total = P,
      width = 80
    )
    max_len <- max(nchar(proteins))
  }
  
  for (i in seq_len(P)) {
    protein <- proteins[i]
    
    if (show_progress) {
      pb$tick(tokens = list(
        what = sprintf(paste0("%-", max_len, "s"), protein)
      ))
    }
    
    chromatogram <- cleaned[protein, ]
    points <- sum(!is.na(profile_matrix[protein, ]))
    
    gaussian <- choose_gaussians_maple(
      chromatogram = chromatogram,
      points = points,
      max_gaussians = max_gaussians,
      criterion = criterion,
      max_iterations = max_iterations,
      min_R_squared = min_R_squared,
      use_R_squared = use_R_squared,
      use_RMSE = use_RMSE,
      max_RMSE = max_RMSE,
      pass_rule = pass_rule,
      optimize_by = optimize_by,
      method = method,
      filter_gaussians_center = filter_gaussians_center,
      filter_gaussians_height = filter_gaussians_height,
      filter_gaussians_variance_min = filter_gaussians_variance_min,
      filter_gaussians_variance_max = filter_gaussians_variance_max
    )
    
    gaussians[[protein]] <- gaussian
  }
  
  return(gaussians)
}


# Assign functions to the PrInCE namespace.
# This is required because several helper functions used above
# are internal to PrInCE.
environment(build_gaussians_maple) <- asNamespace("PrInCE")
ns <- asNamespace("PrInCE")

environment(get_fraction_index) <- ns
environment(fit_gaussians_maple) <- ns
environment(choose_gaussians_maple) <- ns
environment(build_gaussians_maple) <- ns
