#' Predictions from an spmboost model
#'
#' Predicts the additive predictor of an [spmboost()] model, by default
#' including the spatial effect, which is kriged to new locations.
#'
#' @param object A fitted model of class `"spmboost"`.
#' @param newdata Optional `data.frame` (or `sf` object with POINT
#'   geometries, transformed to the CRS of the training data) with the
#'   covariates and -- whenever the spatial effect is needed -- the
#'   coordinate columns used in [spmboost()]. If `NULL`, predictions are for
#'   the training data.
#' @param spatial What to do with the spatial effect: `TRUE` (default)
#'   adds the kriged spatial effect, `FALSE` predicts without it (no
#'   coordinates needed), `"only"` returns the kriged spatial effect alone
#'   (final iteration; not combined with `aggregate = "cumsum"`).
#' @param type Scale of the predictions: `"link"` (default), `"response"`, or
#'   `"class"` (binomial only).
#' @param which Optional subset of base-learners (names or indices). Then,
#'   as in [mboost::predict.mboost()], the contributions of these
#'   base-learners are returned, without intercept and spatial effect.
#' @param aggregate `"sum"` (default): prediction at `mstop`; `"cumsum"`:
#'   matrix with one column per iteration `1, ..., mstop` (used by
#'   [spat_cvrisk()]); `"none"`: per-iteration contributions of the
#'   base-learners, as in mboost (without intercept and spatial effect).
#' @param local Optional integer `k`: krige with the `k` nearest training
#'   locations instead of all of them -- an approximation that keeps
#'   prediction cheap for very large data. Used automatically (with
#'   `k = 100`) when the model was fitted with local fitting (N > 3000).
#' @param ... Further arguments passed to [mboost::predict.mboost()].
#'
#' @details The predictor at iteration \eqn{m} is
#'   \deqn{\hat\eta = f^{[m]} + \hat\beta_0^{[m]} + \hat s^{[m]},}
#'   with \eqn{f^{[m]}} the offset plus all base-learner contributions,
#'   \eqn{\hat\beta_0^{[m]}} the intercept of the spatial model and
#'   \eqn{\hat s^{[m]}} the spatial effect, kriged from the training
#'   locations with that iteration's covariance parameters. If the spatial
#'   effect was excluded in an iteration (see `spat_threshold` in
#'   [spmboost()]), \eqn{\hat s^{[m]} = 0}.
#'
#' @return A numeric vector (or a matrix for `aggregate = "cumsum"` and
#'   `"none"`, as in mboost).
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
#' new <- data.frame(sx = c(0.2, 0.8), sy = c(0.5, 0.5), x1 = 0.5)
#' predict(mod, newdata = new)                     # with spatial effect
#' predict(mod, newdata = new, spatial = FALSE)    # without
#' predict(mod, newdata = new, spatial = "only")   # spatial effect alone
#' dim(predict(mod, newdata = new, aggregate = "cumsum"))
#'
#' @seealso [spmboost()], [extract_spatial()]
#' @export
predict.spmboost <- function(object,
                             newdata = NULL,
                             spatial = TRUE,
                             type = c("link", "response", "class"),
                             which = NULL,
                             aggregate = c("sum", "cumsum", "none"),
                             local = NULL,
                             ...) {
  internalpredict.spmboost(object = object,
                           newdata = newdata,
                           spatial = spatial,
                           type = type,
                           which = which,
                           aggregate = aggregate,
                           local = local,
                           ...)
}


#' Fitted values of an spmboost model
#'
#' In-sample fitted values on the response scale, including the spatial
#' effect, i.e. `predict(object, type = "response")`.
#'
#' @param object A fitted model of class `"spmboost"`.
#' @param spatial Include the spatial effect (`TRUE`, default) or not.
#' @param ... Further arguments passed to [predict.spmboost()].
#'
#' @return Numeric vector of fitted values, one per training observation.
#'
#' @seealso [predict.spmboost()]
#' @export
fitted.spmboost <- function(object, spatial = TRUE, ...) {
  predict.spmboost(object, newdata = NULL, spatial = spatial, type = "response", ...)
}


#' @rdname spat_cv
#' @param object A fold matrix of class `"spat_cv_boost"` from `spat_cv()`.
#' @export
autoplot.spat_cv_boost <- function(object, ...) {
  cv_obj <- attr(object, "cv_object")
  if (is.null(cv_obj)) {
    stop("autoplot() is only available for spatial folds (cv_type = 'spatial', 'block', ",
         "'buffered' or 'loo'); random folds have no spatial geometry to plot.")
  }
  autoplot(cv_obj, ...)
}


#' Methods for cross-validation results
#'
#' `mstop()` returns the CV-optimal number of boosting iterations; `plot()`
#' shows the out-of-sample risk of every fold (grey), their average (black)
#' and the optimal iteration (dashed line).
#'
#' @param object,x An object of class `"spat_cv"` from [spat_cvrisk()].
#' @param main Plot title.
#' @param ... Not used (`mstop()`), or further graphical parameters passed
#'   to [graphics::plot()].
#'
#' @return `mstop()`: the CV-optimal number of iterations. `plot()`: called
#'   for its side effect.
#'
#' @seealso [spat_cvrisk()]
#' @name spat_cv-methods
NULL

#' @rdname spat_cv-methods
#' @export
mstop.spat_cv <- function(object, ...) {
  object$m_opt
}

#' @rdname spat_cv-methods
#' @export
plot.spat_cv <- function(x, main = "CV Risks", ...) {
  avg_risks <- x$avg_risks
  RISKMAT   <- x$cv_risks
  iter      <- length(avg_risks)

  plot(x = seq_len(iter),
       y = avg_risks,
       type = "l",
       ylim = range(RISKMAT),
       main = main,
       xlab = "Number of Iterations",
       ylab = "out-of-bag Risk",
       ...)

  for (i in seq_len(ncol(RISKMAT))) {
    lines(seq_len(iter), y = RISKMAT[, i], col = "gray")
  }
  lines(x = seq_len(iter), y = avg_risks)
  abline(v = x$m_opt, lty = "dashed", col = "gray")
}
