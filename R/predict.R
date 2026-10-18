#' Predict from a locally-weighted PLS model
#'
#' For every row of `new_data`, a weighted PLS model is fitted around that
#' sample and used to predict it.
#'
#' @param object A `lwpls_fit` object created by [lwpls_fit()].
#' @param new_data A data frame or matrix of new predictors.
#' @param type A single character string, or `NULL`. One of:
#'   * `"numeric"` for numeric predictions (regression; the default).
#'   * `"class"` for hard class predictions (classification; the default).
#'   * `"prob"` for class probabilities (classification).
#'   * `"raw"` for an array with the predictions of the local models with
#'     `1, ..., num_comp` components (see Value).
#' @param num_comp The number of PLS components to use. Defaults to the value
#'   used in [lwpls_fit()]. Since the model is lazy, any value can be used;
#'   values larger than the maximum possible number of components give the
#'   same predictions as that maximum.
#' @param ... Not used, but required for extensibility.
#'
#' @return
#' A tibble of predictions with one row per row of `new_data`, following the
#' tidymodels conventions: `.pred` (or `.pred_{outcome}` for multiple
#' outcomes), `.pred_class`, or `.pred_{level}` columns.
#'
#' For `type = "raw"`, a numeric array of dimension
#' `nrow(new_data) x q x num_comp`, where `q` is the number of outcomes (or of
#' classes, for the predicted class indicators), on the original outcome
#' scale.
#'
#' Rows of `new_data` with missing values get missing predictions.
#'
#' @seealso [lwpls_fit()], [multi_predict._lwpls_fit()]
#' @examples
#' train <- mtcars[-(1:5), ]
#' test <- mtcars[1:5, ]
#'
#' fit <- lwpls_fit(mpg ~ ., data = train, num_comp = 3, localization = 0.5)
#' predict(fit, test)
#' predict(fit, test, num_comp = 1)
#' predict(fit, test, type = "raw")
#' @export
predict.lwpls_fit <- function(object,
                              new_data,
                              type = NULL,
                              num_comp = NULL,
                              ...) {
  rlang::check_dots_empty()
  type <- check_pred_type(object, type, allow_raw = TRUE)
  num_comp <- check_pred_num_comp(object, num_comp, single = TRUE)

  forged <- hardhat::forge(new_data, object$blueprint)
  # Robust models are refitted for each number of components: only compute
  # the ones needed.
  comps <- if (type == "raw") seq_len(num_comp) else num_comp
  raw <- lwpls_predict_array(object, forged$predictors, num_comp, comps)
  if (type == "raw") {
    return(raw)
  }

  res <- format_predictions(object, raw_slice(raw, num_comp), type)
  hardhat::validate_prediction_size(res, new_data)
  res
}

#' Predictions for several numbers of components
#'
#' `multi_predict()` returns the predictions of a fitted [lwpls()] model for
#' several values of `num_comp` at once. All values are obtained from a
#' single pass over the data.
#'
#' @param object A [parsnip::model_fit] object created from a [lwpls()]
#'   specification.
#' @param new_data A rectangular data object, such as a data frame.
#' @param type A single character value or `NULL`. Possible values are
#'   `"numeric"`, `"class"`, or `"prob"`. When `NULL`, `"numeric"` is used for
#'   regression and `"class"` for classification.
#' @param num_comp An integer vector with the numbers of components. Defaults
#'   to the value used to fit the model.
#' @param ... Not currently used.
#' @return A tibble with the same number of rows as `new_data` and a list
#'   column `.pred` of tibbles, each with a `num_comp` column and the
#'   prediction column(s).
#' @examples
#' fit <-
#'   lwpls(num_comp = 4) |>
#'   parsnip::set_mode("regression") |>
#'   parsnip::fit(mpg ~ ., data = mtcars[-(1:5), ])
#'
#' preds <- parsnip::multi_predict(fit, mtcars[1:5, ], num_comp = 1:4)
#' preds
#' preds$.pred[[1]]
#' @export
multi_predict._lwpls_fit <- function(object,
                                     new_data,
                                     type = NULL,
                                     num_comp = NULL,
                                     ...) {
  rlang::check_dots_empty()

  fit <- object$fit
  type <- check_pred_type(fit, type, allow_raw = FALSE)
  num_comp <- check_pred_num_comp(fit, num_comp, single = FALSE)

  new_data <- parsnip::prepare_data(object, new_data)
  forged <- hardhat::forge(new_data, fit$blueprint)
  raw <- lwpls_predict_array(fit, forged$predictors, max(num_comp), num_comp)
  n <- dim(raw)[1]

  res <- lapply(num_comp, function(a) {
    pred <- format_predictions(fit, raw_slice(raw, a), type)
    tibble::add_column(pred, num_comp = a, .before = 1)
  })
  res <- vctrs::vec_rbind(!!!res)
  rows <- rep(seq_len(n), times = length(num_comp))
  res <- res[order(rows, res$num_comp), , drop = FALSE]

  tibble::tibble(.pred = vctrs::vec_split(res, sort(rows))$val)
}

# ------------------------------------------------------------------------------

# Predictions for 1, ..., num_comp components as an m x q x num_comp array on
# the original outcome scale. Rows with missing values give NA. For robust
# models, only the numbers of components in `comps` are computed; the other
# slices are NA.
lwpls_predict_array <- function(object, predictors, num_comp, comps = seq_len(num_comp)) {
  x <- as.matrix(predictors[names(object$x_center)])
  storage.mode(x) <- "double"
  m <- nrow(x)
  q <- length(object$y_center)

  res <- array(
    NA_real_,
    dim = c(m, q, num_comp),
    dimnames = list(NULL, object$outcome_names, seq_len(num_comp))
  )
  complete <- rowSums(!is.finite(x)) == 0
  if (!any(complete)) {
    return(res)
  }

  n_comp <- min(num_comp, object$max_comp)
  new_x <- standardize(x[complete, , drop = FALSE], object$x_center, object$x_scale)
  args <- kernel_args(object, new_x)

  if (args$robust > 0 || args$sparsity > 0) {
    pred <- lwpls_general_cpp(
      x = object$x,
      y = object$y,
      new_x = new_x,
      dist_x = args$dist_x,
      dist_new = args$dist_new,
      comps = sort(unique(pmin(comps, n_comp))),
      localization = object$localization,
      neighbors = object$neighbors,
      sparsity = args$sparsity,
      robust = args$robust,
      fair_c = args$fair_c,
      hampel_probs = args$hampel_probs,
      max_iter = args$max_iter,
      classification = object$mode == "classification",
      tol = lwpls_tol
    )
    pred[is.nan(pred)] <- NA_real_
  } else {
    pred <- lwpls_predict_cpp(
      x = object$x,
      y = object$y,
      new_x = new_x,
      dist_x = args$dist_x,
      dist_new = args$dist_new,
      num_comp = n_comp,
      localization = object$localization,
      neighbors = object$neighbors,
      tol = lwpls_tol
    )
  }

  n_slices <- dim(pred)[3]
  for (j in seq_len(q)) {
    res[complete, j, seq_len(n_slices)] <-
      pred[, j, ] * object$y_scale[[j]] + object$y_center[[j]]
  }
  # No component can be added beyond the maximum: repeat the last one.
  if (num_comp > n_comp) {
    res[complete, , seq(n_comp + 1, num_comp)] <- res[complete, , n_comp]
  }
  res
}

# Settings of the C++ kernels for a fitted model and standardized queries.
kernel_args <- function(object, new_x) {
  robust <- isTRUE(object$robust)
  fair <- !robust || object$weight_function == "fair"
  if (is.null(object$projection)) {
    dist_x <- object$x
    dist_new <- new_x
  } else {
    dist_x <- object$x %*% object$projection
    dist_new <- new_x %*% object$projection
  }
  list(
    dist_x = dist_x,
    dist_new = dist_new,
    sparsity = if (is.null(object$sparsity)) 0 else object$sparsity,
    robust = if (!robust) 0L else if (fair) 1L else 2L,
    fair_c = if (robust && fair) object$robust_constant else 4,
    hampel_probs = if (robust && !fair) object$robust_constant else c(0.95, 0.975, 0.999),
    max_iter = if (robust) object$max_iter else 1L
  )
}

# Relative tolerance used to stop extracting components from a local model
# whose predictors or covariances are exhausted.
lwpls_tol <- 1e-12

raw_slice <- function(raw, num_comp) {
  d <- dim(raw)
  matrix(raw[, , num_comp], nrow = d[1], ncol = d[2], dimnames = list(NULL, dimnames(raw)[[2]]))
}

format_predictions <- function(object, pred, type) {
  switch(
    type,
    numeric = {
      if (ncol(pred) == 1) {
        hardhat::spruce_numeric(pred[, 1])
      } else {
        cols <- lapply(seq_len(ncol(pred)), function(j) pred[, j])
        names(cols) <- colnames(pred)
        hardhat::spruce_numeric_multiple(!!!cols)
      }
    },
    class = {
      cls <- object$lvl[max.col(pred, ties.method = "first")]
      hardhat::spruce_class(factor(cls, levels = object$lvl))
    },
    prob = hardhat::spruce_prob(object$lvl, indicators_to_prob(pred))
  )
}

# Predicted class indicators sum to one; truncating them to [0, 1] and
# renormalizing gives proper probabilities.
indicators_to_prob <- function(pred) {
  prob <- pmin(pmax(pred, 0), 1)
  total <- rowSums(prob)
  degenerate <- !is.na(total) & total <= 0
  prob[degenerate, ] <- 1
  total[degenerate] <- ncol(prob)
  prob <- prob / total
  dimnames(prob) <- NULL
  prob
}

check_pred_type <- function(object, type, allow_raw, call = rlang::caller_env()) {
  if (is.null(type)) {
    return(if (object$mode == "regression") "numeric" else "class")
  }
  valid <- if (object$mode == "regression") "numeric" else c("class", "prob")
  if (allow_raw) {
    valid <- c(valid, "raw")
  }
  rlang::arg_match0(type, valid, arg_nm = "type", error_call = call)
}

check_pred_num_comp <- function(object, num_comp, single, call = rlang::caller_env()) {
  if (is.null(num_comp)) {
    return(object$num_comp)
  }
  ok <- is.numeric(num_comp) && length(num_comp) >= 1 &&
    !anyNA(num_comp) && all(num_comp >= 1) && all(num_comp == round(num_comp))
  if (single) {
    ok <- ok && length(num_comp) == 1
  }
  if (!ok) {
    what <- if (single) "a single whole number" else "a vector of whole numbers"
    cli::cli_abort("{.arg num_comp} must be {what} >= 1.", call = call)
  }
  sort(unique(as.integer(num_comp)))
}
