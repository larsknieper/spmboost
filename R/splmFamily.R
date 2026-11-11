as.Family <- function(object, ...)
  UseMethod("as.Family")

## --- Pool of spatially coherent partitions (for spmodel's local fitting) ---
build_local_index_pool <- function(data, xcoord, ycoord, k, pool_size = 20) {
  N    <- nrow(data)
  pool <- matrix(NA_integer_, nrow = N, ncol = pool_size)

  for (b in seq_len(pool_size)) {
    SCVM <- spat_cv(data, k = k, coords = c(xcoord, ycoord))
    # SCVM: N x k, 1 = training, 0 = test
    stopifnot(all(rowSums(SCVM == 0) == 1))  # one 0 per line expected
    pool[, b] <- apply(SCVM, 1, which.min)
  }
  pool
}


## as.Family.splm() builds an mboost Family() boosted spatial models.
## Each boosting iteration re-fits splm() on the offset-adjusted data, then
## caches kriging for prediction.
##
## spat_threshold (t_s): after every spatial fit, the spatial share of the
## residual variance, de / (de + ie) = sigma_s^2 / (sigma_s^2 + sigma^2), is
## computed. If it is below t_s, no spatial effect is included. The spatial model
## is still fitted every iteration, so if the share stays below 
## t_s throughout, the result is a fixed-effects-only model.
as.Family.splm <- function(object,
                           xcoord = NULL,
                           ycoord = NULL,
                           spat_threshold   = 0,     # t_s; 0 = spatial effect always included
                           local_threshold  = 3000,  # activate spmodel's local fitting above this N
                           local_k          = 10,    # no. of spatial groups (~ N/local_k for each group)
                           local_pool_size  = 20,    # how many partitions in advance
                           local_var_adjust = "theoretical",
                           local_parallel   = FALSE,
                           optim_args = list(reltol = 1e-2),
                           cold_restart_every = 10, # every k-th iteration, drop the warm-started
                           verbose          = TRUE,
                           ...) {

  model      <- object
  iter_count <- 0L   # counts ngradient() calls; drives cold_restart_every

  if (inherits(model, "spglm")) {
    fam <- eval(parse(text = model$family))
  } else if (inherits(model, "splm")) {
    fam <- gaussian()
  } else {
    stop("as.Family.splm() expects an object fitted with spmodel::splm() or spmodel::spglm().")
  }

  mf <- model.frame(model)
  N  <- nrow(mf)
  Y  <- model.response(mf)

  ## --- Per-iteration state, updated by ngradient() via `<<-` ---
  ## Defined here, so `<<-` writes into THIS closure (one per fitted model),
  ## never into the user's global environment. mboost stores a copy of the
  ## state after every iteration (nuisance[[m]] <- family@nuisance()), so
  ## mod$nuisance()[[m]] / extract_spatial(mod, m) give each iteration's values,
  ## independently of other models fitted before or after (e.g. CV refits).
  ff            <- 0
  fitted_spat   <- rep(0, N)
  SpatCovPars   <- NULL
  resid_i       <- NULL
  gamma_i       <- NULL
  spatial_share <- NA_real_
  spatial_in    <- TRUE

  use_local <- N > local_threshold

  index_pool <- NULL
  if (use_local) {
    index_pool <- build_local_index_pool(model$obdata, xcoord = xcoord, ycoord = ycoord,
                                         k = local_k, pool_size = local_pool_size)
    if (verbose) {
      message(sprintf(paste0("N = %d > local_threshold = %d: %d precomputed spatial partitions ",
                             "(k = %d groups) were generated; one is drawn at random each iteration."),
                      N, local_threshold, local_pool_size, local_k))
    }
  }

  ## --- Metadata needed later for out-of-sample (by-hand) prediction ---
  coords_train <- data.frame(
    xcoord = model$obdata[[xcoord]],
    ycoord = model$obdata[[ycoord]]
  )
  # Coordinates never change across iterations -- only the covariance PARAMETERS
  # do. Caching the distance matrix here means each iteration only has to
  # recompute the (cheap, elementwise) correlation transform, not the O(n^2)
  # distances themselves.
  h_train <- if (!use_local) as.matrix(dist(coords_train)) else NULL

  if (use_local && verbose) {
    message("Note: with local (big-data) fitting active, exact global dual-kriging weights ",
            "are NOT precomputed per iteration (would defeat the purpose of the local ",
            "approximation). The predict function will fall back to k-nearest-neighbour ",
            "local kriging for the spatial term; pass `local = <k>` to predict() ",
            "to control the neighbourhood size.")
  }

  grad_formula <- update(model$formula, . ~ . + offset(f))

  risk <- function(y, f, w = 1) {
    if (length(f) == 1) f <- rep(f, N)
    # Reconstruct eta from the CURRENT f, not a cached model. 
    lp <- f + ff + fitted_spat
    if (fam$family == "gaussian") {
      return(sum(w * (Y - lp)^2))
    }
    mu  <- fam$linkinv(lp)
    dev <- sum(fam$dev.resids(Y, mu, w))
    fam$aic(y = Y, n = 1, mu = mu, wt = w, dev = dev) / 2
  }

  ngradient <- function(y, f, w = 1) {
    if (length(f) == 1) f <- rep(f, N)
    if (length(w) == 1) w <- rep(w, N)
    stopifnot(NROW(f) == N)

    iter_count <<- iter_count + 1L
    do_cold_restart <- isTRUE(cold_restart_every > 0) &&
      (iter_count %% cold_restart_every == 0)

    if (do_cold_restart) {
      init_spcov <- NULL
      if (verbose) {
        message(sprintf("[splm-family] iteration %d: cold-starting spcov (skipping warm start).",
                        iter_count))
      }
    } else {
      init_spcov <- list(initial = model$coefficients$spcov,
                         is_known = model$is_known$spcov)
      attr(init_spcov, "class") <- attributes(model$coefficients$spcov)$class
    }

    ## --- Offset-adjusted GLS refit + cached dual-kriging weights for prediction.
    data_with_f   <- model$obdata
    data_with_f$f <- as.numeric(f)

    local_arg <- FALSE
    if (use_local) {
      col_pick  <- sample.int(local_pool_size, 1)
      local_arg <- list(index      = index_pool[, col_pick],
                        var_adjust = local_var_adjust,
                        parallel   = local_parallel)
    }

    # optim_args is passed straight through to update() -> splm() -> stats::optim()
    # (splm()'s `...` is documented as passed to esv()/stats::optim()).
    # Any further splm() arguments given via spmboost(spmodel_control = ...)
    # are part of the stored call and therefore carried along by update().
    model <<- withCallingHandlers(
      do.call(update, c(list(object        = model,
                             formula       = grad_formula,
                             data          = data_with_f,
                             spcov_initial = init_spcov,
                             local         = local_arg,
                             evaluate      = TRUE),
                        optim_args)),
      message = function(m) {
        if (grepl("spcov_initial overriding spcov_type", conditionMessage(m), fixed = TRUE)) {
          invokeRestart("muffleMessage")
        }
      })

    # Intercept-only formulas are assumed throughout (no other fixed covariates),
    ff_gls <- coef(model)["(Intercept)"]
    # Raw residual y - f (spmodel's response residuals are marginal.
    r <- model$residuals$response + ff_gls

    SpatCovPars <<- model$coefficients$spcov

    # --- Spatial selection: share of the residual variance that is spatial ---
    de <- SpatCovPars[["de"]]
    ie <- SpatCovPars[["ie"]]
    spatial_share <<- if (isTRUE(de + ie > 0)) as.numeric(de / (de + ie)) else 0
    spatial_in    <<- spatial_share >= spat_threshold

    if (!spatial_in) {
      # No spatial effect in this iteration: zero spatial fit, no kriging
      # weights, and iid errors -> the intercept is the plain mean.
      ff          <<- mean(r)
      fitted_spat <<- rep(0, N)
      resid_i     <<- as.numeric(r - ff)
      gamma_i     <<- NULL
      return(resid_i)                              # y - f - ff
    }

    ff          <<- ff_gls
    fitted_spat <<- model$fitted$spcov$de
    resid_i     <<- model$residuals$response   # y - f - ff (marginal)

    # --- Precompute dual-kriging weights ONCE per iteration ---
    # gamma = Sigma^{-1} %*% resid_i.
    # Skipped when local (big-data) fitting is active, since forming/solving the
    # full n x n Sigma here would defeat the point of the local approximation;
    # the predict function falls back to local (k-NN) kriging in that case.
    gamma_i <<- NULL
    if (!use_local) {
      R_train <- .spcov_corr(h_train, SpatCovPars)
      Sigma_o <- SpatCovPars[["de"]] * R_train
      diag(Sigma_o) <- diag(Sigma_o) + SpatCovPars[["ie"]]
      L       <- tryCatch(chol(Sigma_o), error = function(e) NULL)
      if (!is.null(L)) {
        gamma_i <<- backsolve(L, backsolve(L, resid_i, transpose = TRUE))
      } else {
        # numerically near-singular Sigma_o (e.g. ie ~ 0): fall back
        gamma_i <<- as.numeric(solve(Sigma_o, resid_i))
      }
    }

    return(as.numeric(resid_i - fitted_spat))     # y - f - ff - s_hat
  }

  Fam <- Family(ngradient = ngradient,
                risk = risk,
                check_y = function(y) rep(TRUE, N),
                offset = function(y, w) 0,
                nuisance = function() list(ff            = ff,
                                           SpatCovPars   = SpatCovPars,
                                           fitted_spat   = fitted_spat,
                                           resid         = resid_i,
                                           gamma         = gamma_i,
                                           spatial_share = spatial_share,
                                           spatial_in    = spatial_in,
                                           local_active  = use_local),
                name = "splm-family",
                response = function(f) fam$linkinv(f)
  )

  # One-time metadata for the predict function (see internal_functions.R).
  attr(Fam, "spmodel_meta") <- list(
    coords_train   = coords_train,
    xcoord         = xcoord,
    ycoord         = ycoord,
    spat_threshold = spat_threshold,
    base_family    = fam
  )

  Fam
}
