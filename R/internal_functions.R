## Prediction for spmboost objects.

## Exact correlation R(h) for every spcov_type splm() 
## `theta` is a spcov_params-like object (named vector with
## de, ie, range, and possibly extra)
## `h` is a distance matrix or vector.
.spcov_corr <- function(h, theta) {
  type <- class(theta)[1]
  type <- sub("^spcov_", "", type)   # normalize regardless of exact class prefix
  
  phi <- theta[["range"]]
  eta <- h / phi
  ind <- function(R) { R[h > phi] <- 0; R }   # I{h <= phi} cutoff, used by several types
  
  switch(type,
         exponential    = exp(-eta),
         gaussian       = exp(-eta^2),
         spherical      = ind(1 - 1.5 * eta + 0.5 * eta^3),
         triangular     = ind(1 - eta),
         circular       = {
           m <- pmin(eta, 1)
           ind(1 - (2 / pi) * (m * sqrt(1 - m^2) + asin(m)))
         },
         cubic          = ind(1 - 7 * eta^2 + 8.75 * eta^3 - 3.5 * eta^5 + 0.75 * eta^7),
         pentaspherical = ind(1 - 1.875 * eta + 1.25 * eta^3 - 0.375 * eta^5),
         cosine         = cos(eta),
         wave           = {
           R <- sin(eta) / eta
           R[h == 0] <- 1
           R
         },
         jbessel        = besselJ(h * phi, nu = 0),
         gravity        = (1 + eta^2)^(-1/2),
         rquad          = (1 + eta^2)^(-1),
         magnetic       = (1 + eta^2)^(-3/2),
         matern         = {
           xi    <- theta[["extra"]]
           alpha <- sqrt(2 * xi * eta)
           R <- (2^(1 - xi) / gamma(xi)) * alpha^xi * besselK(alpha, xi)
           R[h == 0] <- 1
           R
         },
         cauchy         = (1 + eta^2)^(-theta[["extra"]]),
         pexponential   = exp(-h^theta[["extra"]] / phi),
         none = , ie    = matrix(0, nrow = NROW(h), ncol = NCOL(h)),
         stop(sprintf(
           paste0("Correlation function for spcov_type '%s' is not implemented in .spcov_corr(). ",
                  "Add it using the R-form table in vignette('technical', 'spmodel')."), type))
  )
}

## Euclidean cross-distance matrix between two coordinate data.frames
.cross_dist <- function(A, B) {
  sqrt(outer(A$xcoord, B$xcoord, "-")^2 + outer(A$ycoord, B$ycoord, "-")^2)
}

## Spatial (kriging) term at new coordinates for ONE boosting iteration.
##  - Fast/exact path: gamma was precomputed at fit time (non-local branch) ->
##    just a cross-covariance times gamma, O(n_train * n_new).
##  - Local/fallback path: k-nearest-neighbour kriging, recomputed per point.
##    Used automatically when gamma is NULL (big-data / local fitting was
##    active), or when the caller passes `local = k` to force it (e.g. to
##    keep prediction on very large new datasets cheap even in the exact
##    branch, at the cost of being an approximation).
.krige_term <- function(coords_new, coords_train, theta, gamma, resid, local = NULL) {
  de <- theta[["de"]]
  
  if (!is.null(gamma) && is.null(local)) {
    R0 <- .spcov_corr(.cross_dist(coords_new, coords_train), theta)
    return(as.numeric((de * R0) %*% gamma))
  }
  
  if (is.null(resid)) {
    stop("Local kriging fallback requires the per-iteration residual vector ",
         "(nuisance()$resid); make sure as.Family.splm() stores it.")
  }
  
  k     <- if (is.null(local)) min(100L, nrow(coords_train)) else local
  n_new <- nrow(coords_new)
  ie    <- theta[["ie"]]
  out   <- numeric(n_new)
  
  for (i in seq_len(n_new)) {
    d_all <- sqrt((coords_train$xcoord - coords_new$xcoord[i])^2 +
                    (coords_train$ycoord - coords_new$ycoord[i])^2)
    idx <- order(d_all)[seq_len(min(k, length(d_all)))]
    
    h_loc     <- as.matrix(dist(coords_train[idx, c("xcoord", "ycoord")]))
    Sigma_loc <- de * .spcov_corr(h_loc, theta) + ie * diag(length(idx))
    c0        <- de * .spcov_corr(d_all[idx], theta)
    
    gamma_loc <- solve(Sigma_loc, resid[idx])
    out[i]    <- sum(c0 * gamma_loc)
  }
  out
}

## Fixed-effects term (Intercept) for one iteration.
.fixed_term <- function(n_new, ff) {
  rep(ff, n_new)
}


## Turn an sf object (POINT geometries) into a plain data.frame that carries
## the coordinates as ordinary columns `xcoord` and `ycoord`. Used for both
## the training data (spmboost()) and newdata (predict()). If `crs` is given
## and differs from the object's CRS, the object is transformed first, so new
## locations are always on the same coordinate scale as the training data.
.sf_to_df <- function(data, xcoord, ycoord, crs = NULL) {
  if (!is.null(crs) && !is.na(crs) && !is.na(st_crs(data)) && st_crs(data) != crs) {
    data <- st_transform(data, crs)
  }
  geom_type <- as.character(st_geometry_type(data))
  if (any(geom_type != "POINT")) {
    warning("Not all geometries are POINTs; using their centroids (sf::st_centroid()).",
            call. = FALSE)
    data <- suppressWarnings(st_centroid(data))
  }
  if (isTRUE(st_is_longlat(data))) {
    warning("data has a geographic (longitude/latitude) CRS: distances are computed as ",
            "Euclidean distances in degrees. Consider projecting first (sf::st_transform()).",
            call. = FALSE)
  }
  coords <- st_coordinates(data)
  out    <- st_drop_geometry(data)
  out[[xcoord]] <- unname(coords[, 1])
  out[[ycoord]] <- unname(coords[, 2])
  list(data = out, crs = st_crs(data))
}


internalpredict.spmboost <- function(object, newdata = NULL, spatial = TRUE,
                                     type = c("link", "response", "class"),
                                     which = NULL, aggregate = c("sum", "cumsum", "none"),
                                     local = NULL, ...) {
  
  type      <- match.arg(type)
  aggregate <- match.arg(aggregate)
  
  meta <- attr(object$family, "spmodel_meta")
  if (is.null(meta)) {
    stop("No spmodel metadata found on object$family. Refit using as.Family.splm() ",
         "from splmFamily.R.")
  }
  
  if (!is.null(newdata) && inherits(newdata, "sf")) {
    newdata <- .sf_to_df(newdata, meta$xcoord, meta$ycoord, crs = object$crs)$data
  }
  
  ## single base-learners / per-iteration contributions: no spatial or intercept term
  if (!identical(spatial, "only") && (aggregate == "none" || !is.null(which))) {
    return(predict.mboost(object = object, newdata = newdata, type = type,
                          which = which, aggregate = aggregate, ...))
  }
  
  coords_train <- meta$coords_train
  need_coords  <- identical(spatial, "only") || isTRUE(spatial)
  coords_from  <- function(nd) {
    if (is.null(nd[[meta$xcoord]]) || is.null(nd[[meta$ycoord]])) {
      stop(sprintf(paste0("newdata is missing the coordinate column(s) '%s' and/or '%s' ",
                          "(needed for the spatial effect; use spatial = FALSE to predict without it)."),
                   meta$xcoord, meta$ycoord))
    }
    data.frame(xcoord = nd[[meta$xcoord]], ycoord = nd[[meta$ycoord]])
  }
  coords_new <- NULL
  if (need_coords) {
    coords_new <- if (is.null(newdata)) coords_train else coords_from(newdata)
  }
  n_new <- if (is.null(newdata)) nrow(coords_train) else nrow(newdata)
  
  # Spatial term of one iteration
  spatial_term <- function(nu) {
    if (isFALSE(nu$spatial_in)) return(numeric(n_new))
    .krige_term(coords_new, coords_train, nu$SpatCovPars, nu$gamma, nu$resid, local)
  }
  
  ## spatial effect only
  if (identical(spatial, "only")) {
    if (aggregate == "cumsum") stop("aggregate = 'cumsum' is not supported for spatial = 'only'.")
    return(spatial_term(object$nuisance()[[object$mstop()]]))
  }
  
  ## boosted (non-spatial) part, via mboost's own predict method
  preds <- predict.mboost(object = object, newdata = newdata, type = "link",
                          which = NULL, aggregate = aggregate, ...)
  
  add_iteration <- function(i) {
    nu  <- object$nuisance()[[i]]
    val <- .fixed_term(n_new, nu$ff)
    if (isTRUE(spatial)) val <- val + spatial_term(nu)
    val
  }
  
  if (aggregate == "sum") {
    preds <- preds + add_iteration(object$mstop())
  } else if (aggregate == "cumsum") {
    if (is.list(preds)) stop("Unexpected list-type prediction for aggregate = 'cumsum'.")
    for (i in seq_len(object$mstop())) preds[, i] <- preds[, i] + add_iteration(i)
  }
  
  if (type %in% c("response", "class")) {
    fam   <- meta$base_family
    preds <- fam$linkinv(preds)
    if (type == "class") {
      if (fam$family == "binomial") {
        preds <- round(preds)
      } else {
        stop("Class-type predictions are only implemented for binomial family.")
      }
    }
  }
  preds
}