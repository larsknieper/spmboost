#' Model-based boosting with a Gaussian Process spatial effect
#'
#' Fits a componentwise gradient boosting model (via [mboost::mboost()]) in
#' which a Gaussian-process spatial effect is re-estimated in every boosting
#' iteration by [spmodel::splm()].
#'
#' @param boost_formula Formula for the boosted (non-spatial) part, as in
#'   [mboost::mboost()].
#' @param data A `data.frame` with coordinate columns, or an `sf` object with
#'   POINT geometries (other geometries are replaced by their centroids). An
#'   `sf` object is converted to a `data.frame` once, up front: its coordinates
#'   become ordinary columns named `xcoord`/`ycoord` (default `".xcoord"`,
#'   `".ycoord"`), and its CRS is stored so that `sf` newdata in `predict()`
#'   is transformed to the same CRS. Distances are Euclidean, so use a
#'   projected CRS.
#' @param spcov_type Spatial covariance type passed to [spmodel::splm()],
#'   e.g. `"gaussian"` (default), `"exponential"`, `"matern"`.
#' @param xcoord,ycoord Names of the x- and y-coordinate columns in `data`.
#'   Required for a `data.frame`; optional for `sf` data (then they name the
#'   coordinate columns that are created from the geometry).
#' @param estmethod Estimation method for the spatial covariance parameters,
#'   passed to [spmodel::splm()]: one of `"reml"` (default), `"ml"`,
#'   `"sv-wls"`, `"sv-cl"`.
#' @param spat_threshold Threshold \eqn{t_s \in [0, 1)} on the spatial share of
#'   the residual variance, \eqn{\hat\sigma_s^2 / (\hat\sigma_s^2 + \hat\sigma^2)}
#'   (partial sill / (partial sill + nugget)). In every boosting iteration in
#'   which the share falls below `spat_threshold`, no spatial effect is
#'   included: it is set to zero in the negative gradient, the risk and the
#'   predictor, so the base-learners are fitted to the residual without it, and
#'   the intercept is the plain mean. The spatial model is still fitted in
#'   every iteration; if the share stays below `spat_threshold` throughout, the
#'   result is a fixed-effects-only model. The default `0` always includes the
#'   spatial effect.
#' @param spmodel_control Named list of further arguments for
#'   [spmodel::splm()]: `spcov_initial`, `weights` (splm's semivariogram
#'   weights, e.g. `"cressie"` -- not to be confused with the boosting
#'   `weights` below), `anisotropy`, `random`, `randcov_initial`,
#'   `partition_factor`, `local`, `range_constrain`, `ddf`, or arguments
#'   passed on to [stats::optim()] (e.g. `method`, `control`). Merged over the
#'   package default `list(range_constrain = TRUE)`; set an element to `NULL`
#'   to drop that default and use splm()'s own. `formula`, `data`,
#'   `spcov_type`, `xcoord`, `ycoord` and `estmethod` are set by spmboost()
#'   itself and must not appear here. Note that splm() ignores `spcov_type`
#'   when `spcov_initial` is given, and that `spcov_initial` and `local` only
#'   affect the initial fit -- from the first boosting iteration on, both are
#'   managed internally (warm starts, and local fitting for large N).
#' @param na.action,weights,offset,control,oobweights,baselearner,...
#'   Arguments passed to [mboost::mboost()]. If `offset` is `NULL`, the link
#'   of the mean response is used.
#' @param family Response distribution. Currently only `gaussian()` with
#'   identity link is supported (as `family` object, function or name);
#'   other families give an error.
#' @param verbose Logical; if `TRUE`, print progress messages of the spatial
#'   fits (cold restarts of the covariance optimisation, local fitting for
#'   large N). Default `FALSE`.
#'
#' @return An object of class `c("spmboost", "mboost")`. Besides the usual
#'   mboost methods, see [predict.spmboost()] for predictions including the
#'   kriged spatial effect and [extract_spatial()] for the spatial fit.
#'
#' @examples
#' set.seed(1)
#' n   <- 100
#' dat <- data.frame(sx = runif(n), sy = runif(n),
#'                   x1 = runif(n), x2 = rnorm(n))
#' # spatial effect: Gaussian process with exponential covariance
#' D <- as.matrix(dist(dat[, c("sx", "sy")]))
#' s <- drop(t(chol(exp(-D / 0.3))) %*% rnorm(n))
#' dat$y <- 1 + sin(2 * pi * dat$x1) + 0.5 * dat$x2 + s + rnorm(n, sd = 0.3)
#'
#' mod <- spmboost(y ~ bbs(x1) + bols(x2), data = dat,
#'                 spcov_type = "exponential",
#'                 xcoord = "sx", ycoord = "sy",
#'                 control = boost_control(mstop = 30))
#'
#' # estimated covariance parameters (de = partial sill, ie = nugget) and
#' # spatial share of the residual variance
#' extract_spatial(mod)[c("SpatCov", "share")]
#'
#' # selected base-learners and partial effect of x1
#' table(selected(mod))
#' plot(mod, which = "x1")
#'
#' # predictions at new locations: full predictor and spatial effect only
#' new <- data.frame(sx = c(0.2, 0.8), sy = c(0.5, 0.5), x1 = 0.5, x2 = 0)
#' predict(mod, newdata = new)
#' predict(mod, newdata = new, spatial = "only")
#'
#' @seealso [spat_cv()] and [spat_cvrisk()] for choosing `mstop`,
#'   [extract_spatial()]
#' @export
spmboost <- function(boost_formula,
                     data = list(),
                     spcov_type = "gaussian",
                     xcoord = NULL,
                     ycoord = NULL,
                     estmethod = "reml",
                     spat_threshold = 0,
                     spmodel_control = list(),
                     na.action = na.omit,
                     weights = NULL,
                     offset = NULL,
                     family = gaussian(),
                     control = boost_control(),
                     oobweights = NULL,
                     baselearner = c("bbs", "bols", "btree", "bss", "bns"),
                     verbose = FALSE,
                     ...) {
  
  if (inherits(family, "boost_family")) {
    stop("Please use a glm-family.")
  }
  
  family <- switch(
    class(family)[1],
    "character" = get(family, mode = "function")(),
    "function"  = family(),
    "family"    = family,
    stop("Invalid family argument")
  )

  # The spatial model is fitted with spmodel::splm() (Gaussian, identity link)
  # and the negative gradient is the Gaussian one -- other families would run
  # without error but silently fit squared-error boosting on the link scale.
  if (family$family != "gaussian" || family$link != "identity") {
    stop("spmboost() currently supports only gaussian() with identity link; got ",
         family$family, "(link = \"", family$link, "\").", call. = FALSE)
  }
  
  ## --- coordinates / sf -----------------------------------------------------
  crs <- NULL
  if (inherits(data, "sf")) {
    if (is.null(xcoord)) xcoord <- ".xcoord"
    if (is.null(ycoord)) ycoord <- ".ycoord"
    sf_prep <- .sf_to_df(data, xcoord, ycoord)
    data    <- sf_prep$data
    crs     <- sf_prep$crs
  } else if (is.null(xcoord) || is.null(ycoord)) {
    stop("Both `xcoord` and `ycoord` (names of the coordinate columns in `data`) must be supplied ",
         "unless `data` is an sf object.")
  }
  
  estmethod <- match.arg(estmethod, c("reml", "ml", "sv-wls", "sv-cl"))
  
  if (!is.numeric(spat_threshold) || length(spat_threshold) != 1 ||
      spat_threshold < 0 || spat_threshold >= 1) {
    stop("`spat_threshold` must be a single number in [0, 1).")
  }
  
  ## --- further splm() arguments --------------------------------------------
  if (!is.list(spmodel_control)) {
    stop("`spmodel_control` must be a (named) list.")
  }
  if (length(spmodel_control) > 0 &&
      (is.null(names(spmodel_control)) || any(names(spmodel_control) == ""))) {
    stop("All elements of `spmodel_control` must be named.")
  }
  reserved <- c("formula", "data", "spcov_type", "xcoord", "ycoord", "estmethod")
  clash    <- intersect(names(spmodel_control), reserved)
  if (length(clash) > 0) {
    stop("The following must be given as arguments of spmboost() itself, not inside ",
         "`spmodel_control`: ", paste(clash, collapse = ", "))
  }
  spmodel_control <- modifyList(list(range_constrain = TRUE), spmodel_control)
  
  # Extract the left-hand side of the complete formula
  lhs <- all.vars(boost_formula)[1]
  
  if (is.null(offset)) {
    offset <- family$linkfun(mean(data[[lhs]]))
  }
  
  # Intercept-only formula for the spatial model
  spm_formula <- as.formula(paste(lhs, "~ 1"))
  
  internal_spm <- do.call(
    splm,
    c(list(formula    = spm_formula,
           data       = data,
           spcov_type = spcov_type,
           xcoord     = xcoord,
           ycoord     = ycoord,
           estmethod  = estmethod),
      spmodel_control)
  )
  
  spm_family <- as.Family.splm(internal_spm,
                               xcoord         = xcoord,
                               ycoord         = ycoord,
                               spat_threshold = spat_threshold,
                               verbose        = verbose)
  
  res <- mboost(formula = boost_formula,
                data = data,
                na.action = na.action,
                weights = weights,
                offset = offset,
                family = spm_family,
                control = control,
                oobweights = oobweights,
                baselearner = baselearner,
                ...)
  class(res) <- c("spmboost", class(res))
  
  res$formula         <- boost_formula
  res$spm_formula     <- spm_formula
  res$sp_family       <- function() family
  res$spcov_type      <- spcov_type
  res$xcoord          <- xcoord
  res$ycoord          <- ycoord
  res$crs             <- crs        # NULL unless data was sf
  res$estmethod       <- estmethod
  res$spat_threshold  <- spat_threshold
  res$spmodel_control <- spmodel_control
  res$def_bl          <- baselearner
  res$data            <- data       # always a data.frame (sf already converted)
  res$control         <- control
  res$verbose         <- verbose
  return(res)
}


#' Extract the estimated spatial effect of an spmboost model
#'
#' Returns the spatial part of a fitted [spmboost()] model for one boosting
#' iteration: the fitted spatial effect at the training locations, the
#' estimated spatial covariance parameters, and whether the spatial effect was
#' included in that iteration (see `spat_threshold` in [spmboost()]).
#'
#' @param object A fitted model of class `"spmboost"`.
#' @param iteration Boosting iteration to extract, between 1 and
#'   `mstop(object)`. Defaults to the final iteration, `mstop(object)`.
#'
#' @return A list with elements
#'   \describe{
#'     \item{`fitted`}{Numeric vector (one value per training observation) of
#'       the fitted spatial effect \eqn{\hat s} at the training locations. All
#'       zeros if the spatial effect was not included in this iteration.}
#'     \item{`SpatCov`}{Estimated spatial covariance parameters from
#'       [spmodel::splm()] (class of the covariance type): `de` (partial sill,
#'       \eqn{\hat\sigma_s^2}), `ie` (nugget, \eqn{\hat\sigma^2}), `range`, plus
#'       `rotate`/`scale` (anisotropy) and `extra` where applicable. Both `de`
#'       and `ie` are variances. Reported even if the spatial effect was not
#'       included.}
#'     \item{`share`}{Spatial share of the residual variance,
#'       \eqn{\hat\sigma_s^2 / (\hat\sigma_s^2 + \hat\sigma^2)}.}
#'     \item{`included`}{Logical: was the spatial effect included in this
#'       iteration, i.e. `share >= spat_threshold`?}
#'   }
#'
#' @details For the spatial effect at new locations use
#'   `predict(object, newdata, spatial = "only")`, which krigs it.
#'
#' @examples
#' set.seed(1)
#' n   <- 100
#' dat <- data.frame(sx = runif(n), sy = runif(n), x1 = runif(n))
#' D   <- as.matrix(dist(dat[, c("sx", "sy")]))
#' dat$y <- sin(2 * pi * dat$x1) + drop(t(chol(exp(-D / 0.3))) %*% rnorm(n)) +
#'   rnorm(n, sd = 0.3)
#' mod <- spmboost(y ~ bbs(x1), data = dat, spcov_type = "exponential",
#'                 xcoord = "sx", ycoord = "sy",
#'                 control = boost_control(mstop = 20))
#'
#' sp <- extract_spatial(mod)
#' sp$SpatCov       # partial sill (de), nugget (ie), range
#' sp$share         # spatial share of the residual variance
#'
#' # estimated partial sill along the boosting path
#' plot(sapply(1:20, function(m) extract_spatial(mod, m)$SpatCov[["de"]]),
#'      type = "l", xlab = "iteration", ylab = "partial sill")
#'
#' @seealso [spmboost()], [predict.spmboost()]
#' @export
extract_spatial <- function(object, iteration = mstop(object)) {
  nu  <- object$nuisance()[[iteration]]
  ret <- list()
  ret$fitted   <- nu$fitted_spat
  ret$SpatCov  <- nu$SpatCovPars
  ret$share    <- nu$spatial_share   # de / (de + ie)
  ret$included <- nu$spatial_in      # share >= spat_threshold
  return(ret)
}


#' Create cross-validation folds for spmboost models
#'
#' Builds a fold matrix for [spat_cvrisk()]: spatial folds via
#' \pkg{spatialsample}, or ordinary random k-fold CV via [mboost::cv()].
#'
#' @param data A `data.frame` (then `coords` is required for the spatial
#'   types) or an `sf` object.
#' @param k Number of folds.
#' @param coords Character vector of length 2: names of the x- and
#'   y-coordinate columns in `data`. Ignored for `sf` data and for
#'   `cv_type = "random"`.
#' @param cv_type Type of cross-validation:
#'   \describe{
#'     \item{`"spatial"`}{Spatial clustering CV: locations are clustered into
#'       `k` groups ([spatialsample::spatial_clustering_cv()]).}
#'     \item{`"block"`}{Spatial block CV: the domain is cut into blocks that
#'       are assigned to folds ([spatialsample::spatial_block_cv()]).}
#'     \item{`"buffered"`}{Random k-fold CV with an exclusion buffer around
#'       each test fold ([spatialsample::spatial_buffer_vfold_cv()]); see
#'       `buffer`.}
#'     \item{`"loo"`}{Leave-location-out CV: groups given by `loo_group` are
#'       held out together ([spatialsample::spatial_leave_location_out_cv()]).}
#'     \item{`"random"`}{Ordinary (non-spatial) random k-fold CV, folds drawn
#'       by `mboost::cv(type = "kfold")`. Needs no coordinates.}
#'   }
#' @param crs Coordinate reference system used when converting a
#'   `data.frame` to `sf` (ignored for `sf` data and `cv_type = "random"`).
#'   The default `NA` uses plain Euclidean distances between the coordinates,
#'   matching the distances used by [spmboost()]. If the coordinates are
#'   longitude/latitude, set `crs = 4326` (or project them first).
#' @param loo_group Name of the grouping column for `cv_type = "loo"`.
#' @param radius,buffer For `cv_type = "buffered"` (both required there,
#'   in the units of the coordinates): observations within `radius` of a
#'   randomly chosen test observation are added to its test fold, and
#'   observations within `buffer` of the test fold are excluded from
#'   training. Excluded observations are coded `-1` in the returned matrix.
#' @param ... Further arguments passed to the \pkg{spatialsample} function
#'   (used for `"buffered"` and `"loo"`).
#'
#' @return A numeric `N x k` matrix of class `"spat_cv_boost"`, one column per
#'   fold, with `1` = training and `0` = test observation (and `-1` =
#'   excluded, for the buffer of `cv_type = "buffered"`). Every observation
#'   is in exactly one test fold, except for `"buffered"` with `radius > 0`,
#'   where test folds can overlap. For the spatial types, the underlying
#'   \pkg{spatialsample} object is attached as attribute `"cv_object"`, which
#'   `autoplot()` uses to plot the folds (not available for `"random"`).
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(sx = runif(200), sy = runif(200), y = rnorm(200))
#'
#' # spatial clustering CV and ordinary random CV, 5 folds each
#' f_sp <- spat_cv(dat, k = 5, coords = c("sx", "sy"))
#' f_rd <- spat_cv(dat, k = 5, cv_type = "random")
#' colSums(f_sp == 0)   # test-fold sizes
#'
#' # buffered CV: points within 0.1 of a test fold are left out of training
#' f_bf <- spat_cv(dat, k = 5, coords = c("sx", "sy"), cv_type = "buffered",
#'                 radius = 0.05, buffer = 0.1)
#' table(f_bf[, 1])     # -1 = excluded, 0 = test, 1 = training
#'
#' \donttest{
#' # plot the spatial folds
#' spatialsample::autoplot(f_sp)
#' }
#'
#' @seealso [spat_cvrisk()]
#' @export
spat_cv <- function(data,
                    k = 10,
                    coords = NULL,
                    cv_type = c("spatial", "block", "buffered", "loo", "random"),
                    crs = NA,
                    loo_group = NULL,
                    radius = NULL,
                    buffer = NULL,
                    ...
) {
  
  cv_type <- match.arg(cv_type)
  
  if (cv_type == "random") {
    # Ordinary (non-spatial) k-fold CV, folds drawn by mboost's own routine.
    folds_M <- unclass(cv(rep(1, nrow(data)), type = "kfold", B = k))
    attributes(folds_M) <- list(dim = dim(folds_M))
    class(folds_M) <- "spat_cv_boost"
    return(folds_M)
  }
  
  if (!inherits(data, "sf")) {
    if (is.null(coords)) stop("No coordinate-columns specified.")
    data_sf <- st_as_sf(
      data,
      coords = coords,
      crs = crs
    )
  } else {
    data_sf <- data
  }
  
  # crs = NA (the default) is deliberate: Euclidean distances on the given
  # coordinates. spatialsample warns about the missing CRS -- muffle exactly
  # that warning in this case; everything else passes through.
  no_crs <- !inherits(data, "sf") && is.na(crs)
  folds_obj <- withCallingHandlers(switch(
    cv_type,
    spatial  = spatial_clustering_cv(data_sf, v = k),
    block    = spatial_block_cv(data_sf, v = k),
    buffered = spatial_buffer_vfold_cv(data_sf, v = k, radius = radius,
                                       buffer = buffer, ...),
    loo      = spatial_leave_location_out_cv(data_sf, v = k, group = loo_group, ...)
  ),
  warning = function(w) {
    if (no_crs && grepl("coordinate reference system", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
  
  # folds_M: N x k, 1 = training, 0 = test, -1 = excluded (buffer only)
  N <- nrow(data)
  folds_M <- matrix(-1, nrow = N, ncol = k)
  for (i in seq_len(k)) {
    split  <- folds_obj$splits[[i]]
    in_id  <- split$in_id
    out_id <- split$out_id
    if (length(out_id) == 1 && is.na(out_id)) out_id <- setdiff(seq_len(N), in_id)
    folds_M[in_id, i]  <- 1
    folds_M[out_id, i] <- 0
  }
  
  attr(folds_M, "cv_object") <- folds_obj
  class(folds_M) <- "spat_cv_boost"
  return(folds_M)
}


#' Cross-validated risk for choosing mstop of an spmboost model
#'
#' For every fold, refits [spmboost()] on the training part with the settings
#' of `object`, predicts the test part for every iteration `1, ..., mstop(object)`,
#' and evaluates the out-of-sample risk. The iteration with the lowest average
#' risk is the CV-optimal `mstop`.
#'
#' @param object A fitted model of class `"spmboost"`. Its current `mstop`
#'   (also after `object[m]`) is the largest number of iterations evaluated.
#' @param folds Fold matrix, typically from [spat_cv()]: `N x k`, `1` =
#'   training, `0` = test (any other value, e.g. `-1` for buffered CV, =
#'   excluded from both), with rows matching `object$data`.
#' @param cores Number of cores; folds are run in parallel via
#'   [parallel::mclapply()] (forking -- on Windows only `cores = 1` works).
#' @param include_spatial Logical. If `TRUE`, the test predictions include the
#'   kriged spatial effect (offset + base-learners + intercept + \eqn{\hat s}),
#'   which is appropriate when test locations lie within the observed domain
#'   (e.g. random folds). If `FALSE` (default), only the non-spatial part is
#'   used (offset + base-learners + intercept).
#'
#' @details Each refit uses the same `boost_formula`, `spcov_type`,
#'   coordinates, `estmethod`, `spat_threshold`, `spmodel_control`, `offset`,
#'   `family`, `control` and `baselearner` as `object`. Arguments that were
#'   passed to [spmboost()] via `...` or `weights` are not stored and
#'   therefore not used in the refits.
#'
#'   The risk per iteration is the squared error for the Gaussian family and
#'   half the AIC (i.e. the negative log-likelihood) otherwise, divided by the
#'   size of the test fold, so folds of different sizes are comparable.
#'
#'   [mboost::cvrisk()] cannot be used for spmboost models.
#'
#' @return An object of class `"spat_cv"` with elements
#'   \describe{
#'     \item{`cv_risks`}{`mstop x k` matrix of out-of-sample risks per
#'       iteration and fold.}
#'     \item{`avg_risks`}{Average risk per iteration (row means of
#'       `cv_risks`).}
#'     \item{`m_opt`}{CV-optimal number of iterations, also returned by
#'       `mstop()`.}
#'     \item{`folds`}{The fold matrix used.}
#'   }
#'   `plot()` shows the risk paths of all folds and their average.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n   <- 60
#' dat <- data.frame(sx = runif(n), sy = runif(n), x1 = runif(n))
#' D   <- as.matrix(dist(dat[, c("sx", "sy")]))
#' dat$y <- sin(2 * pi * dat$x1) + drop(t(chol(exp(-D / 0.3))) %*% rnorm(n)) +
#'   rnorm(n, sd = 0.3)
#'
#' mod   <- spmboost(y ~ bbs(x1), data = dat, spcov_type = "exponential",
#'                   xcoord = "sx", ycoord = "sy",
#'                   control = boost_control(mstop = 20))
#' folds <- spat_cv(dat, k = 3, cv_type = "random")
#' cvr   <- spat_cvrisk(mod, folds, include_spatial = TRUE)
#' mstop(cvr)
#' plot(cvr)
#'
#' # set the model to the CV-optimal number of iterations
#' # (same result as refitting with mstop = mstop(cvr))
#' mod <- mod[mstop(cvr)]
#' }
#'
#' @seealso [spat_cv()], [spmboost()]
#' @export
spat_cvrisk <- function(object, folds,
                        cores = 1, include_spatial = FALSE) {
  k <- ncol(folds)
  iter <- object$mstop()
  folds_l <- lapply(seq_len(k), function(i) folds[, i])

  # refits run exactly as many iterations as are evaluated -- the current
  # mstop(object), which may differ from the stored control after object[m]
  ctrl <- object$control
  ctrl$mstop <- iter
  
  Y <- object$data[[all.vars(object$formula[[2]])]]
  
  # internal cross-validation function
  cross_val <- function(giv_weights) {
    train_ind <- giv_weights == 1
    test_ind  <- giv_weights == 0
    
    train <- object$data[train_ind, ]
    test  <- object$data[test_ind, ]
    
    train_mod <- spmboost(boost_formula   = object$formula,
                          data            = train,
                          spcov_type      = object$spcov_type,
                          xcoord          = object$xcoord,
                          ycoord          = object$ycoord,
                          estmethod       = object$estmethod,
                          spat_threshold  = object$spat_threshold,
                          spmodel_control = object$spmodel_control,
                          offset          = object$offset,
                          family          = object$sp_family,
                          control         = ctrl,
                          baselearner     = object$def_bl,
                          verbose         = isTRUE(object$verbose))
    
    eta <- predict(train_mod,
                   newdata = test,
                   type = "link",
                   aggregate = "cumsum",
                   spatial = include_spatial)
    
    aic <- numeric(iter)
    
    fam <- object$sp_family()
    
    for (j in seq_len(iter)) {
      mu <- fam$linkinv(eta[, j])
      if (fam$family == "gaussian") {
        aic[j] <- sum((Y[test_ind] - mu)^2) * 2
      } else {
        dev <- sum(fam$dev.resids(y = Y[test_ind], mu = mu, wt = 1))
        aic[j] <- fam$aic(y = Y[test_ind], n = 1, mu = mu, wt = 1, dev = dev)
      }
    }
    return(aic / 2)
  }
  aics <- mclapply(folds_l, cross_val, mc.cores = cores)
  RISKMAT <- matrix(NA, nrow = iter, ncol = k)
  for (i in seq_along(aics)) {
    RISKMAT[, i] <- aics[[i]]
  }
  
  test_sizes <- colSums(unclass(folds) == 0)
  
  RISKMAT <- sweep(RISKMAT, 2, test_sizes, "/")
  
  avg_risks <- rowMeans(RISKMAT)
  m_opt <- which.min(avg_risks)
  
  res <- list(cv_risks = RISKMAT,
              avg_risks = avg_risks,
              m_opt = m_opt,
              folds = folds)
  
  class(res) <- "spat_cv"
  
  return(res)
}