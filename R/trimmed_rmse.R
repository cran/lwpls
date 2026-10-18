#' Trimmed root mean squared error
#'
#' The root mean squared error of the `1 - trim` fraction of the observations
#' with the smallest absolute errors. When the reference values contain
#' outliers, they also appear in the assessment sets of a resampling scheme,
#' where their errors would dominate the usual RMSE and favor non-robust
#' models. The trimmed RMSE leaves them out, like the trimmed standard error
#' of prediction used to tune partial robust M-regression (Serneels et al.,
#' 2005), so it is the metric to use when tuning robust LW-PLS models
#' ([lwpls()] with `robust = TRUE`) on contaminated data.
#'
#' @param data A data frame containing the columns given by `truth` and
#'   `estimate`.
#' @param truth The column identifier for the true results (numeric).
#' @param estimate The column identifier for the predicted results (numeric).
#' @param trim The fraction of observations with the largest absolute errors
#'   to leave out, in `[0, 1)`. The default, 0.2, leaves out 20%.
#' @param na_rm A logical: should missing values be removed?
#' @param case_weights The optional column identifier for case weights.
#' @param ... Not currently used.
#' @return A tibble with columns `.metric`, `.estimator` and `.estimate` and
#'   one row (or one row per group), or, for `trimmed_rmse_vec()`, a single
#'   number.
#' @references
#' Serneels, S., Croux, C., Filzmoser, P. and Van Espen, P. J. (2005). Partial
#' robust M-regression. *Chemometrics and Intelligent Laboratory Systems*,
#' 79(1--2), 55--64. \doi{10.1016/j.chemolab.2005.04.007}
#' @examples
#' df <- data.frame(truth = c(1, 2, 3, 4, 5), estimate = c(1.1, 2.1, 2.9, 4.2, 15))
#' trimmed_rmse(df, truth, estimate)
#' yardstick::rmse(df, truth, estimate)
#'
#' # In a metric set, for example with tune::tune_grid()
#' metrics <- yardstick::metric_set(trimmed_rmse, yardstick::rmse)
#' metrics(df, truth, estimate)
#' @export
trimmed_rmse <- function(data, ...) {
  UseMethod("trimmed_rmse")
}
trimmed_rmse <- yardstick::new_numeric_metric(trimmed_rmse, direction = "minimize")

#' @rdname trimmed_rmse
#' @export
trimmed_rmse.data.frame <- function(data,
                                    truth,
                                    estimate,
                                    trim = 0.2,
                                    na_rm = TRUE,
                                    case_weights = NULL,
                                    ...) {
  yardstick::numeric_metric_summarizer(
    name = "trimmed_rmse",
    fn = trimmed_rmse_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    na_rm = na_rm,
    case_weights = !!rlang::enquo(case_weights),
    fn_options = list(trim = trim)
  )
}

#' @rdname trimmed_rmse
#' @export
trimmed_rmse_vec <- function(truth,
                             estimate,
                             trim = 0.2,
                             na_rm = TRUE,
                             case_weights = NULL,
                             ...) {
  ok <- is.numeric(trim) && length(trim) == 1 && !is.na(trim) && trim >= 0 && trim < 1
  if (!ok) {
    cli::cli_abort("{.arg trim} must be a single number in [0, 1), not {describe(trim)}.")
  }
  yardstick::check_numeric_metric(truth, estimate, case_weights)

  if (na_rm) {
    result <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- result$truth
    estimate <- result$estimate
    case_weights <- result$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }

  trimmed_rmse_impl(truth, estimate, trim, case_weights)
}

# Keeps the observations with the smallest absolute errors that hold a
# fraction 1 - trim of the total weight (at least one observation).
trimmed_rmse_impl <- function(truth, estimate, trim, case_weights) {
  sq_error <- (truth - estimate)^2
  weights <- if (is.null(case_weights)) rep(1, length(sq_error)) else as.numeric(case_weights)
  if (length(sq_error) == 0) {
    return(NA_real_)
  }
  o <- order(sq_error)
  sq_error <- sq_error[o]
  weights <- weights[o]
  keep <- cumsum(weights) <= (1 - trim) * sum(weights) * (1 + 1e-12)
  keep[1] <- TRUE
  sqrt(sum(weights[keep] * sq_error[keep]) / sum(weights[keep]))
}
