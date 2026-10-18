#' Fit a locally-weighted partial least squares model
#'
#' @description
#' `lwpls_fit()` is the computational engine behind [lwpls()]. It can also be
#' used on its own through the matrix, data frame, formula or recipe
#' interfaces.
#'
#' LW-PLS is a lazy learner: fitting validates and standardizes the training
#' data and stores it together with the hyperparameters. The local models are
#' built at prediction time (see [predict.lwpls_fit()]).
#'
#' Sparse (`sparsity > 0`) and robust (`robust = TRUE`) local models are
#' described in [lwpls()].
#'
#' @param x Depending on the context:
#'   * A __data frame__ or __matrix__ of numeric predictors.
#'   * A __recipe__ specifying a set of preprocessing steps created from
#'     [recipes::recipe()].
#' @param y When `x` is a __data frame__ or __matrix__, `y` is the outcome:
#'   * A __factor__ (or a one-column data frame with a factor) for
#'     classification.
#'   * A __numeric vector__, or a __matrix__/__data frame__ with one or more
#'     numeric columns, for regression.
#' @param data When a __recipe__ or __formula__ is used, `data` is a data
#'   frame containing both the predictors and the outcome(s).
#' @param formula A formula specifying the outcome terms on the left-hand
#'   side and the predictor terms on the right-hand side. Multiple numeric
#'   outcomes can be given as `y1 + y2 ~ .`.
#' @param num_comp The number of PLS components in each local model. It is
#'   capped at the maximum possible value, `min(ncol(x), nrow(x) - 1)`.
#' @param localization A positive number: the localization parameter of the
#'   similarity weights (see [lwpls()]). Smaller values give more local
#'   models.
#' @param neighbors Either `NULL` (default), to use all training samples, or
#'   the number of nearest training samples used to build each local model.
#' @param similarity The similarity index used to weight the training
#'   samples: `"euclidean"` (default) or `"covariance"` (covariance-based
#'   LW-PLS, CbLW-PLS). See [lwpls()] for details.
#' @param sparsity A number in `[0, 1)`: the sparsity threshold \eqn{\eta} of
#'   the local models (Hoffmann et al., 2015). Zero (default) gives ordinary
#'   local PLS models.
#' @param scale A logical: should the predictors (and numeric outcomes) be
#'   standardized to unit variance? Predictors are always centered, which
#'   does not affect the predictions. For robust models, the center and scale
#'   are the median and the median absolute deviation (MAD).
#' @param robust A logical: should the local models be robust to outliers
#'   (partial robust M-regression, Serneels et al., 2005)?
#' @param weight_function The weight function of the robust models:
#'   `"fair"` (default) or `"hampel"`.
#' @param robust_constant The tuning constant(s) of the weight function: for
#'   `"fair"`, the constant \eqn{c} of the Fair function (default 4); for
#'   `"hampel"`, three increasing probabilities that define the cutoffs
#'   (default `c(0.95, 0.975, 0.999)`). `NULL` uses the defaults.
#' @param max_iter The maximum number of iterations of the robust models.
#' @param ... Not currently used, but required for extensibility.
#'
#' @return A `lwpls_fit` object, which stores the standardized training data,
#'   the standardization constants and the hyperparameters.
#'
#' @seealso [predict.lwpls_fit()], [lwpls()]
#' @examples
#' train <- mtcars[-(1:5), ]
#' test <- mtcars[1:5, ]
#'
#' # XY interface
#' fit <- lwpls_fit(train[, -1], train$mpg, num_comp = 3, localization = 0.5)
#' fit
#' predict(fit, test[, -1])
#'
#' # Formula interface, with two outcomes
#' fit2 <- lwpls_fit(mpg + qsec ~ ., data = train, num_comp = 3)
#' predict(fit2, test)
#'
#' # Classification
#' fit3 <- lwpls_fit(Species ~ ., data = iris, num_comp = 2)
#' predict(fit3, iris[c(1, 51, 101), ], type = "prob")
#'
#' # Robust and sparse local models
#' fit4 <- lwpls_fit(mpg ~ ., data = train, num_comp = 3, robust = TRUE, sparsity = 0.3)
#' fit4
#' predict(fit4, test)
#' @export
lwpls_fit <- function(x, ...) {
  UseMethod("lwpls_fit")
}

#' @export
#' @rdname lwpls_fit
lwpls_fit.default <- function(x, ...) {
  cli::cli_abort("{.fn lwpls_fit} is not defined for {.obj_type_friendly {x}}.")
}

#' @export
#' @rdname lwpls_fit
lwpls_fit.data.frame <- function(x,
                                 y,
                                 num_comp = 2L,
                                 localization = 1,
                                 neighbors = NULL,
                                 similarity = "euclidean",
                                 sparsity = 0,
                                 scale = TRUE,
                                 robust = FALSE,
                                 weight_function = "fair",
                                 robust_constant = NULL,
                                 max_iter = 30L,
                                 ...) {
  rlang::check_dots_empty()
  check_character_outcome(y)
  processed <- hardhat::mold(x, y)
  lwpls_bridge(
    processed,
    num_comp = num_comp,
    localization = localization,
    neighbors = neighbors,
    similarity = similarity,
    sparsity = sparsity,
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = max_iter
  )
}

#' @export
#' @rdname lwpls_fit
lwpls_fit.matrix <- function(x,
                             y,
                             num_comp = 2L,
                             localization = 1,
                             neighbors = NULL,
                             similarity = "euclidean",
                             sparsity = 0,
                             scale = TRUE,
                             robust = FALSE,
                             weight_function = "fair",
                             robust_constant = NULL,
                             max_iter = 30L,
                             ...) {
  rlang::check_dots_empty()
  check_character_outcome(y)
  processed <- hardhat::mold(x, y)
  lwpls_bridge(
    processed,
    num_comp = num_comp,
    localization = localization,
    neighbors = neighbors,
    similarity = similarity,
    sparsity = sparsity,
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = max_iter
  )
}

#' @export
#' @rdname lwpls_fit
lwpls_fit.formula <- function(formula,
                              data,
                              num_comp = 2L,
                              localization = 1,
                              neighbors = NULL,
                              similarity = "euclidean",
                              sparsity = 0,
                              scale = TRUE,
                              robust = FALSE,
                              weight_function = "fair",
                              robust_constant = NULL,
                              max_iter = 30L,
                              ...) {
  rlang::check_dots_empty()
  # With an intercept, factors get the same reference-cell coding as with
  # parsnip's formula interface. The intercept column is dropped later.
  processed <- hardhat::mold(
    formula,
    data,
    blueprint = hardhat::default_formula_blueprint(
      intercept = TRUE,
      indicators = "traditional"
    )
  )
  lwpls_bridge(
    processed,
    num_comp = num_comp,
    localization = localization,
    neighbors = neighbors,
    similarity = similarity,
    sparsity = sparsity,
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = max_iter
  )
}

#' @export
#' @rdname lwpls_fit
lwpls_fit.recipe <- function(x,
                             data,
                             num_comp = 2L,
                             localization = 1,
                             neighbors = NULL,
                             similarity = "euclidean",
                             sparsity = 0,
                             scale = TRUE,
                             robust = FALSE,
                             weight_function = "fair",
                             robust_constant = NULL,
                             max_iter = 30L,
                             ...) {
  rlang::check_dots_empty()
  processed <- hardhat::mold(x, data)
  lwpls_bridge(
    processed,
    num_comp = num_comp,
    localization = localization,
    neighbors = neighbors,
    similarity = similarity,
    sparsity = sparsity,
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = max_iter
  )
}

# ------------------------------------------------------------------------------
# Bridge: validate the molded data and hyperparameters

lwpls_bridge <- function(processed,
                         num_comp,
                         localization,
                         neighbors,
                         similarity,
                         sparsity,
                         scale,
                         robust,
                         weight_function,
                         robust_constant,
                         max_iter,
                         call = rlang::caller_env()) {
  check_whole(num_comp, "num_comp", min = 1, call = call)
  check_positive(localization, "localization", call = call)
  check_whole(neighbors, "neighbors", min = 1, allow_null = TRUE, call = call)
  check_similarity(similarity, call = call)
  check_sparsity(sparsity, call = call)
  check_bool(scale, "scale", call = call)
  check_bool(robust, "robust", call = call)
  weight_function <- check_weight_function(weight_function, call = call)
  robust_constant <- check_robust_constant(robust_constant, weight_function, call = call)
  check_whole(max_iter, "max_iter", min = 1, call = call)

  predictors <- processed$predictors
  predictors <- predictors[names(predictors) != "(Intercept)"]
  outcomes <- processed$outcomes

  nms <- names(predictors)
  if (any(is.na(nms) | nms == "") || anyDuplicated(nms) > 0) {
    cli::cli_abort(
      "The predictors must have unique, non-empty column names.",
      call = call
    )
  }

  is_numeric <- vapply(predictors, is.numeric, logical(1))
  if (!all(is_numeric)) {
    bad <- names(predictors)[!is_numeric]
    cli::cli_abort(
      c(
        "All predictors must be numeric.",
        "x" = "Non-numeric predictor{?s}: {.var {bad}}.",
        "i" = "Use the formula interface or a recipe to create dummy variables."
      ),
      call = call
    )
  }
  x <- as.matrix(predictors)
  storage.mode(x) <- "double"

  if (nrow(x) < 2) {
    cli::cli_abort("At least two training samples are required.", call = call)
  }
  if (ncol(x) < 1) {
    cli::cli_abort("At least one predictor is required.", call = call)
  }
  if (!all(is.finite(x))) {
    cli::cli_abort(
      "The predictors contain missing or infinite values.",
      call = call
    )
  }

  if (ncol(outcomes) == 1 && is.factor(outcomes[[1]])) {
    mode <- "classification"
    y_raw <- outcomes[[1]]
    if (anyNA(y_raw)) {
      cli::cli_abort("The outcome contains missing values.", call = call)
    }
    lvl <- levels(y_raw)
    if (length(lvl) < 2) {
      cli::cli_abort("The outcome factor must have at least two levels.", call = call)
    }
    y <- class_indicators(y_raw)
    outcome_names <- lvl
  } else if (all(vapply(outcomes, is.numeric, logical(1)))) {
    mode <- "regression"
    lvl <- NULL
    y <- as.matrix(outcomes)
    storage.mode(y) <- "double"
    if (!all(is.finite(y))) {
      cli::cli_abort(
        "The outcome contains missing or infinite values.",
        call = call
      )
    }
    outcome_names <- colnames(y)
  } else {
    cli::cli_abort(
      c(
        "Unsupported outcome type.",
        "i" = "Use a single factor (classification) or one or more numeric
               columns (regression)."
      ),
      call = call
    )
  }

  lwpls_fit_impl(
    x = x,
    y = y,
    mode = mode,
    lvl = lvl,
    outcome_names = outcome_names,
    num_comp = num_comp,
    localization = localization,
    neighbors = neighbors,
    similarity = similarity,
    sparsity = sparsity,
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = max_iter,
    blueprint = processed$blueprint
  )
}

# ------------------------------------------------------------------------------
# Implementation

lwpls_fit_impl <- function(x,
                           y,
                           mode,
                           lvl,
                           outcome_names,
                           num_comp,
                           localization,
                           neighbors,
                           similarity,
                           sparsity,
                           scale,
                           robust,
                           weight_function,
                           robust_constant,
                           max_iter,
                           blueprint) {
  n <- nrow(x)
  p <- ncol(x)

  max_comp <- min(p, n - 1L)
  if (num_comp > max_comp) {
    cli::cli_warn(c(
      "!" = "{.arg num_comp} = {num_comp} is larger than the maximum number of
             components ({max_comp}).",
      "i" = "{max_comp} component{?s} will be used."
    ))
    num_comp <- max_comp
  }

  if (is.null(neighbors)) {
    neighbors <- n
  } else if (neighbors > n) {
    cli::cli_warn(c(
      "!" = "{.arg neighbors} = {neighbors} is larger than the number of
             training samples ({n}).",
      "i" = "All {n} samples will be used."
    ))
    neighbors <- n
  }

  # Robust models are standardized with the median and the MAD, so that
  # outliers do not distort the distances.
  center_fun <- if (robust) column_median else colMeans
  scale_fun <- if (robust) column_mad else column_sd

  x_center <- center_fun(x)
  x_scale <- if (scale) scale_fun(x) else rep(1, p)
  names(x_center) <- names(x_scale) <- colnames(x)

  # Class indicators are only centered by their means: their predictions then
  # sum to one.
  if (mode == "regression") {
    y_center <- center_fun(y)
    y_scale <- if (scale) scale_fun(y) else rep(1, ncol(y))
  } else {
    y_center <- colMeans(y)
    y_scale <- rep(1, ncol(y))
  }

  x <- standardize(x, x_center, x_scale)
  y <- standardize(y, y_center, y_scale)

  # CbLW-PLS measures distances after projecting the samples on the
  # covariance direction(s) X'Y (Hazama & Kano, 2015).
  projection <- NULL
  if (similarity == "covariance") {
    projection <- if (robust) {
      robust_covariance_projection(x, y, mode == "classification")
    } else {
      covariance_projection(x, y)
    }
  }

  hardhat::new_model(
    x = x,
    y = y,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    y_scale = y_scale,
    mode = mode,
    lvl = lvl,
    outcome_names = outcome_names,
    num_comp = as.integer(num_comp),
    max_comp = as.integer(max_comp),
    localization = as.numeric(localization),
    neighbors = as.integer(neighbors),
    similarity = similarity,
    projection = projection,
    sparsity = as.numeric(sparsity),
    scale = scale,
    robust = robust,
    weight_function = weight_function,
    robust_constant = robust_constant,
    max_iter = as.integer(max_iter),
    blueprint = blueprint,
    class = "lwpls_fit"
  )
}

#' @export
print.lwpls_fit <- function(x, ...) {
  n <- nrow(x$x)
  cat("Locally-weighted PLS (", x$mode, ")\n\n", sep = "")
  cat("Training samples:", n, "\n")
  cat("Predictors:      ", ncol(x$x), if (x$scale) "(standardized)", "\n")
  if (x$mode == "classification") {
    cat("Classes:         ", paste(x$lvl, collapse = ", "), "\n")
  } else if (!identical(x$outcome_names, ".outcome")) {
    cat("Outcomes:        ", paste(x$outcome_names, collapse = ", "), "\n")
  }
  cat("Components:      ", x$num_comp, "\n")
  cat("Localization:    ", format(x$localization, digits = 4), "\n")
  cat(
    "Similarity:      ",
    if (x$similarity == "covariance") "covariance-based (CbLW-PLS)" else "Euclidean",
    "\n"
  )
  cat(
    "Neighbors:       ",
    if (x$neighbors >= n) paste0("all (", n, ")") else x$neighbors,
    "\n"
  )
  if (x$sparsity > 0) {
    cat("Sparsity:        ", format(x$sparsity, digits = 4), "(SNIPLS)\n")
  }
  if (x$robust) {
    constant <- if (x$weight_function == "fair") {
      paste0("c = ", format(x$robust_constant, digits = 4))
    } else {
      paste0("p = ", paste(format(x$robust_constant, digits = 4), collapse = ", "))
    }
    fun <- if (x$weight_function == "fair") "Fair" else "Hampel"
    cat("Robust:           ", "PRM, ", fun, " weights (", constant, ")\n", sep = "")
  }
  invisible(x)
}

# ------------------------------------------------------------------------------
# Helpers

check_character_outcome <- function(y, call = rlang::caller_env()) {
  is_chr <- if (is.data.frame(y)) any(vapply(y, is.character, logical(1))) else is.character(y)
  if (is_chr) {
    cli::cli_abort(
      c(
        "Character outcomes are not supported.",
        "i" = "Convert the outcome to a factor for classification."
      ),
      call = call
    )
  }
  invisible(y)
}

# Gamma = X'Y / ||X'Y||, the first PLS weight vector for a single outcome.
# A covariance below sqrt(eps) of its Cauchy-Schwarz bound ||X|| ||Y|| is
# rounding noise: Gamma is then zero, so all samples are equally similar.
covariance_projection <- function(x, y) {
  xy <- crossprod(x, y)
  norm <- sqrt(sum(xy^2))
  if (norm > sqrt(.Machine$double.eps) * sqrt(sum(x^2) * sum(y^2))) {
    xy / norm
  } else {
    xy * 0
  }
}

# Robust version of Gamma: X'WY / ||X'WY||, with the initial weights of partial
# robust M-regression (Fair function of the distances to the center and, for
# numeric outcomes, of the absolute centered outcomes).
robust_covariance_projection <- function(x, y, classification) {
  w <- fair_weight(row_norms(x) / positive_median(row_norms(x)))
  if (!classification) {
    yc <- sweep(y, 2, column_median(y))
    w <- w * fair_weight(row_norms(yc) / positive_median(row_norms(yc)))
  }
  covariance_projection(x * sqrt(w), y * sqrt(w))
}

fair_weight <- function(z, c = 4) {
  1 / (1 + abs(z / c))^2
}

row_norms <- function(x) {
  sqrt(rowSums(x^2))
}

# Median of non-negative values used as a scale, ignoring zeros if needed.
positive_median <- function(x) {
  res <- stats::median(x)
  if (res > 0) {
    return(res)
  }
  if (any(x > 0)) stats::median(x[x > 0]) else 1
}

column_median <- function(x) {
  apply(x, 2, stats::median)
}

column_mad <- function(x) {
  res <- apply(x, 2, stats::mad)
  bad <- !is.finite(res) | res <= .Machine$double.eps * pmax(abs(column_median(x)), 1)
  res[bad] <- column_sd(x[, bad, drop = FALSE])
  res
}

class_indicators <- function(y) {
  lvl <- levels(y)
  res <- matrix(0, nrow = length(y), ncol = length(lvl), dimnames = list(NULL, lvl))
  res[cbind(seq_along(y), as.integer(y))] <- 1
  res
}

column_sd <- function(x) {
  res <- apply(x, 2, stats::sd)
  res[!is.finite(res) | res <= .Machine$double.eps * pmax(abs(colMeans(x)), 1)] <- 1
  res
}

standardize <- function(x, center, scale) {
  x <- sweep(x, 2, center, check.margin = FALSE)
  x <- sweep(x, 2, scale, "/", check.margin = FALSE)
  dimnames(x) <- NULL
  x
}
