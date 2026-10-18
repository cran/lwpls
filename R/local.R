#' Local models of LW-PLS predictions
#'
#' @description
#' `lwpls_local()` returns the local models that [lwpls()] builds to predict
#' new samples, with the diagnostics that matter in practice, such as for
#' spectroscopic calibration:
#'
#' * **Reliability** of each prediction: the effective number of neighbors,
#'   the distance to the nearest training sample, and the Hotelling
#'   \eqn{T^2} and \eqn{Q} statistics of the query in its local model,
#'   compared with theoretical and empirical limits and ranked among the
#'   training samples of the local model. Queries beyond the limits are
#'   extrapolations of their local model. See [plot_reliability()].
#' * **Regression coefficients** and variable importance in projection (VIP)
#'   of each local model, on the original scale of the predictors and
#'   outcomes. See [plot_coefficients()].
#' * **Selected predictors** of sparse local models. See [plot_selection()].
#' * **Robust weights** of the training samples, averaged over the local
#'   models in which they take part. See [plot_robust_weights()].
#'
#' @param object A fitted [lwpls()] model: a `lwpls_fit` object from
#'   [lwpls_fit()], a parsnip model fit, or a fitted workflow.
#' @param new_data A data frame or matrix of new samples (the queries).
#' @param num_comp The number of components of the local models. Defaults to
#'   the value used to fit the model.
#' @param wavelength An optional numeric vector with the position of each
#'   predictor on the spectral axis (wavelengths, wavenumbers, ...). By
#'   default, the positions are the numbers at the end of the predictor names
#'   (for example `nm_1650` or `x_001`) when they are all present and unique,
#'   and the column numbers otherwise.
#' @param limits The limits that define the `status` and the ratios of the
#'   reliability table: `"empirical"` (default) or `"theoretical"`. Both are
#'   always computed; see Details.
#' @param level The coverage of the limits, and the percentile rank beyond
#'   which a query is flagged. Defaults to 0.95.
#' @param ... Not currently used.
#'
#' @details
#' For a query with local model weights \eqn{w_i} (the similarity weights,
#' multiplied by the robust weights for robust models), the effective number
#' of neighbors is \eqn{(\sum w_i)^2 / \sum w_i^2}{(sum w)^2 / sum w^2}.
#'
#' The \eqn{T^2} statistic of a sample is its squared distance to the
#' weighted mean of the local scores, scaled by their weighted variances, and
#' its \eqn{Q} statistic is the squared norm of its residual (on the
#' standardized scale). They are computed for the query and for the training
#' samples of its local model, and the query is compared with them in three
#' ways:
#'
#' * **Theoretical limits** (`t2_limit_theoretical`, `q_limit_theoretical`):
#'   the \eqn{T^2} limit of a new observation for a model with \eqn{a}
#'   components fitted on \eqn{n} samples,
#'   \eqn{a (n^2 - 1) / (n (n - a)) F_{level}(a, n - a)}{a (n^2 - 1) / (n (n - a)) F(level; a, n - a)},
#'   where \eqn{n} is the effective number of samples of the local model, and
#'   the Box approximation of the weighted distribution of the \eqn{Q} of the
#'   training samples (Nomikos and MacGregor, 1995). They assume normally
#'   distributed scores and residuals, and they are derived for ordinary,
#'   unweighted models.
#' * **Empirical limits** (`t2_limit_empirical`, `q_limit_empirical`): the
#'   `level` quantiles of the \eqn{T^2} and \eqn{Q} of the training samples
#'   of the local model, weighted by their weights in the model. They make no
#'   distributional assumption. The training samples were used to fit the
#'   model, so their statistics are slightly smaller than those of new
#'   samples, which makes these limits slightly conservative (more queries
#'   are flagged).
#' * **Percentile ranks** (`t2_percentile`, `q_percentile`): the weighted
#'   fraction of the training samples of the local model whose statistic is
#'   not larger than that of the query. A rank above `level` is equivalent to
#'   exceeding the empirical limit, on a continuous scale.
#'
#' A query is flagged as `"few neighbors"` when its local model has too few
#' effective samples to assess it (no more than \eqn{a + 1}), whatever the
#' limits.
#'
#' @return An object of class `lwpls_local`, a list with:
#' * `reliability`: a tibble with one row per query: `.row`, the
#'   prediction(s), `n_eff`, `nearest_distance`, `num_comp`, `t2`,
#'   `t2_limit_theoretical`, `t2_limit_empirical`, `t2_percentile`, `q`,
#'   `q_limit_theoretical`, `q_limit_empirical`, `q_percentile`, `t2_ratio`
#'   and `q_ratio` (the statistics divided by the limits chosen by `limits`),
#'   and `status`.
#' * `coefficients`: a tibble with one row per query, outcome and predictor:
#'   `.row`, `outcome`, `predictor`, `position`, `coefficient`, `vip` and
#'   `selected`.
#' * `selection`: a tibble with one row per component and predictor:
#'   `component`, `predictor`, `position`, `frequency` (the fraction of the
#'   local models with this component whose weight vector selects the
#'   predictor) and `queries` (the number of those local models). Sparse
#'   models select predictors per component; a predictor is used by a local
#'   model (`selected` in `coefficients`) if any component selects it.
#' * `training`: a tibble with one row per training sample: `.row`, the
#'   outcome(s), `influence` (the sum of its similarity weights over the
#'   queries) and `robust_weight` (the average of its robust weights,
#'   weighted by the similarity weights).
#' * `info`: the settings of the model and of the spectral axis.
#'
#' @references
#' Nomikos, P. and MacGregor, J. F. (1995). Multivariate SPC charts for
#' monitoring batch processes. *Technometrics*, 37(1), 41--59.
#' \doi{10.1080/00401706.1995.10485888}
#'
#' @seealso [plot_reliability()], [plot_coefficients()],
#'   [plot_robust_weights()], [plot_selection()]
#' @examples
#' fit <- lwpls_fit(mpg ~ ., data = mtcars[-(1:5), ], num_comp = 3)
#' local <- lwpls_local(fit, mtcars[1:5, ])
#' local
#' local$reliability
#'
#' # Theoretical limits at 99%
#' lwpls_local(fit, mtcars[1:5, ], limits = "theoretical", level = 0.99)
#' @export
lwpls_local <- function(object, ...) {
  UseMethod("lwpls_local")
}

#' @export
#' @rdname lwpls_local
lwpls_local.default <- function(object, ...) {
  cli::cli_abort(
    "{.fn lwpls_local} is not defined for {.obj_type_friendly {object}}."
  )
}

#' @export
#' @rdname lwpls_local
lwpls_local.lwpls_fit <- function(object,
                                  new_data,
                                  num_comp = NULL,
                                  wavelength = NULL,
                                  limits = c("empirical", "theoretical"),
                                  level = 0.95,
                                  ...) {
  rlang::check_dots_empty()
  num_comp <- check_pred_num_comp(object, num_comp, single = TRUE)
  limits <- rlang::arg_match(limits)
  check_level(level)
  forged <- hardhat::forge(new_data, object$blueprint)
  local_models(object, forged$predictors, num_comp, wavelength, limits, level)
}

#' @export
#' @rdname lwpls_local
lwpls_local.model_fit <- function(object, new_data, ...) {
  if (!inherits(object$fit, "lwpls_fit")) {
    cli::cli_abort("{.arg object} must be a fitted {.fn lwpls} model.")
  }
  lwpls_local(object$fit, parsnip::prepare_data(object, new_data), ...)
}

#' @export
#' @rdname lwpls_local
lwpls_local.workflow <- function(object, new_data, ...) {
  rlang::check_installed("workflows")
  mold <- workflows::extract_mold(object)
  forged <- hardhat::forge(new_data, mold$blueprint)
  lwpls_local(workflows::extract_fit_parsnip(object), forged$predictors, ...)
}

# ------------------------------------------------------------------------------

local_models <- function(object, predictors, num_comp, wavelength, limits, level) {
  predictor_names <- names(object$x_center)
  axis <- spectral_axis(predictor_names, wavelength)

  x <- as.matrix(predictors[predictor_names])
  storage.mode(x) <- "double"
  m <- nrow(x)
  p <- ncol(x)
  q <- length(object$y_center)
  n <- nrow(object$x)
  complete <- rowSums(!is.finite(x)) == 0
  n_comp <- min(num_comp, object$max_comp)

  pred <- matrix(NA_real_, m, q, dimnames = list(NULL, object$outcome_names))
  coef <- array(NA_real_, c(m, p, q))
  vip <- matrix(NA_real_, m, p)
  selected <- matrix(NA, m, p)
  rel <- matrix(NA_real_, m, 11)
  omega_sum <- rep(0, n)
  omega_robust_sum <- rep(0, n)
  selected_comp <- matrix(0, p, n_comp)
  comp_count <- rep(0, n_comp)

  if (any(complete)) {
    new_x <- standardize(x[complete, , drop = FALSE], object$x_center, object$x_scale)
    args <- kernel_args(object, new_x)
    res <- lwpls_local_cpp(
      x = object$x,
      y = object$y,
      new_x = new_x,
      dist_x = args$dist_x,
      dist_new = args$dist_new,
      num_comp = n_comp,
      localization = object$localization,
      neighbors = object$neighbors,
      sparsity = args$sparsity,
      robust = args$robust,
      fair_c = args$fair_c,
      hampel_probs = args$hampel_probs,
      max_iter = args$max_iter,
      classification = object$mode == "classification",
      tol = lwpls_tol,
      level = level
    )
    pred[complete, ] <- sweep(sweep(res$pred, 2, object$y_scale, "*"), 2, object$y_center, "+")
    # Coefficients on the original scales of the predictors and outcomes
    for (j in seq_len(q)) {
      b <- matrix(res$coef[, , j], nrow = sum(complete), ncol = p)
      coef[complete, , j] <- sweep(b, 2, object$x_scale, "/") * object$y_scale[[j]]
    }
    vip[complete, ] <- res$vip
    selected[complete, ] <- res$selected == 1
    rel[complete, ] <- res$reliability
    # Quantities that cannot be computed are NaN in the C++ code
    vip[is.nan(vip)] <- NA_real_
    rel[is.nan(rel)] <- NA_real_
    omega_sum <- as.vector(res$omega_sum)
    omega_robust_sum <- as.vector(res$omega_robust_sum)
    selected_comp <- res$selected_comp
    comp_count <- as.vector(res$comp_count)
  }

  selection <- tibble::tibble(
    component = rep(seq_len(n_comp), each = p),
    predictor = factor(rep(names(object$x_center), times = n_comp), levels = names(object$x_center)),
    position = rep(axis$position, times = n_comp),
    frequency = as.vector(sweep(selected_comp, 2, pmax(comp_count, 1), "/")),
    queries = rep(as.integer(comp_count), each = p)
  )
  selection <- selection[selection$queries > 0, , drop = FALSE]

  new_local(object, pred, coef, vip, selected, rel, omega_sum, omega_robust_sum,
            selection, axis, n_comp, limits, level)
}

new_local <- function(object, pred, coef, vip, selected, rel, omega_sum,
                      omega_robust_sum, selection, axis, n_comp, limits, level) {
  m <- nrow(pred)
  p <- dim(coef)[2]
  q <- dim(coef)[3]
  predictor_names <- names(object$x_center)

  predictions <- if (object$mode == "classification") {
    format_predictions(object, pred, "class")
  } else {
    format_predictions(object, pred, "numeric")
  }

  # Too few effective samples to assess the query: no theoretical T2 limit
  few <- !is.na(rel[, 1]) & is.na(rel[, 5])
  reliability <- tibble::tibble(
    .row = seq_len(m),
    predictions,
    n_eff = rel[, 1],
    nearest_distance = rel[, 2],
    num_comp = as.integer(rel[, 3]),
    t2 = rel[, 4],
    t2_limit_theoretical = rel[, 5],
    t2_limit_empirical = ifelse(few, NA_real_, rel[, 8]),
    t2_percentile = ifelse(few, NA_real_, rel[, 10]),
    q = rel[, 6],
    q_limit_theoretical = rel[, 7],
    q_limit_empirical = ifelse(few, NA_real_, rel[, 9]),
    q_percentile = ifelse(few, NA_real_, rel[, 11])
  )
  ratios <- reliability_ratios(reliability, limits)
  reliability$t2_ratio <- ratios$t2
  reliability$q_ratio <- ratios$q
  reliability$status <- reliability_status(reliability, ratios)

  outcome_names <- object$outcome_names
  coefficients <- tibble::tibble(
    .row = rep(rep(seq_len(m), times = p), times = q),
    outcome = factor(rep(outcome_names, each = m * p), levels = outcome_names),
    predictor = factor(rep(rep(predictor_names, each = m), times = q), levels = predictor_names),
    position = rep(rep(axis$position, each = m), times = q),
    coefficient = as.vector(coef),
    vip = rep(as.vector(vip), times = q),
    selected = rep(as.vector(selected), times = q)
  )
  coefficients <- coefficients[order(coefficients$.row, coefficients$outcome), , drop = FALSE]

  y <- sweep(sweep(object$y, 2, object$y_scale, "*"), 2, object$y_center, "+")
  outcomes <- if (object$mode == "classification") {
    tibble::tibble(.class = factor(object$lvl[max.col(y, ties.method = "first")], levels = object$lvl))
  } else {
    tibble::as_tibble(as.data.frame(y, optional = TRUE), .name_repair = "minimal") |>
      stats::setNames(object$outcome_names)
  }
  training <- tibble::tibble(
    .row = seq_len(nrow(object$x)),
    outcomes,
    influence = omega_sum,
    robust_weight = ifelse(omega_sum > 0, omega_robust_sum / omega_sum, NA_real_)
  )

  structure(
    list(
      reliability = reliability,
      coefficients = coefficients,
      selection = selection,
      training = training,
      info = list(
        num_comp = n_comp,
        limits = limits,
        level = level,
        mode = object$mode,
        outcome_names = outcome_names,
        robust = isTRUE(object$robust),
        sparsity = if (is.null(object$sparsity)) 0 else object$sparsity,
        axis_label = axis$label,
        position = axis$position,
        center = unname(object$x_center)
      )
    ),
    class = "lwpls_local"
  )
}

#' @export
print.lwpls_local <- function(x, ...) {
  rel <- x$reliability
  cat("LW-PLS local models\n\n")
  cat("Queries:          ", nrow(rel), "\n")
  cat("Components:       ", x$info$num_comp, "\n")
  counts <- table(rel$status)
  counts <- counts[counts > 0]
  if (length(counts) > 0) {
    cat(
      "Reliability:      ", paste0(counts, " ", names(counts), collapse = ", "),
      paste0("(", x$info$limits, " limits, ", format(100 * x$info$level), "%)"), "\n"
    )
  }
  if (x$info$sparsity > 0) {
    n_sel <- tapply(x$coefficients$selected, x$coefficients$.row, function(s) sum(s, na.rm = TRUE)) /
      length(x$info$outcome_names)
    cat(
      "Selected:         ", stats::median(n_sel, na.rm = TRUE), "of",
      length(x$info$position), "predictors (median)\n"
    )
  }
  if (x$info$robust) {
    rw <- x$training[!is.na(x$training$robust_weight), , drop = FALSE]
    lowest <- rw[order(rw$robust_weight)[seq_len(min(3, nrow(rw)))], , drop = FALSE]
    cat(
      "Robust weights:   ", "median", format(stats::median(rw$robust_weight), digits = 2),
      "| lowest:", paste0("#", lowest$.row, " (", format(lowest$robust_weight, digits = 2), ")", collapse = ", "),
      "\n"
    )
  }
  invisible(x)
}

# Statistics divided by the limits of the given kind ("theoretical" or
# "empirical"); a zero statistic gives a zero ratio.
reliability_ratios <- function(reliability, limits) {
  ratio <- function(stat, limit) ifelse(stat == 0, 0, stat / limit)
  list(
    t2 = ratio(reliability$t2, reliability[[paste0("t2_limit_", limits)]]),
    q = ratio(reliability$q, reliability[[paste0("q_limit_", limits)]])
  )
}

reliability_status <- function(reliability, ratios) {
  high_t2 <- !is.na(ratios$t2) & ratios$t2 > 1
  high_q <- !is.na(ratios$q) & ratios$q > 1
  status <- ifelse(
    high_t2 & high_q,
    "high T2 and Q",
    ifelse(high_t2, "high T2", ifelse(high_q, "high Q", "inside"))
  )
  status[is.na(reliability$t2_limit_theoretical)] <- "few neighbors"
  status[is.na(reliability$n_eff)] <- NA
  factor(status, levels = c("inside", "high T2", "high Q", "high T2 and Q", "few neighbors"))
}

check_level <- function(x, arg = "level", call = rlang::caller_env()) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0 && x < 1
  if (!ok) {
    cli::cli_abort("{.arg {arg}} must be a single number in (0, 1), not {describe(x)}.", call = call)
  }
  invisible(x)
}

# Positions of the predictors on the spectral axis.
spectral_axis <- function(predictors, wavelength = NULL, call = rlang::caller_env()) {
  p <- length(predictors)
  if (!is.null(wavelength)) {
    ok <- is.numeric(wavelength) && length(wavelength) == p &&
      all(is.finite(wavelength)) && !anyDuplicated(wavelength)
    if (!ok) {
      cli::cli_abort(
        "{.arg wavelength} must be {p} unique finite number{?s}, one per predictor.",
        call = call
      )
    }
    return(list(position = as.numeric(wavelength), label = "Wavelength"))
  }
  number <- regmatches(
    predictors,
    regexpr("[0-9]+(\\.[0-9]+)?(?=[^0-9]*$)", predictors, perl = TRUE)
  )
  if (length(number) == p) {
    position <- as.numeric(number)
    if (!anyDuplicated(position)) {
      label <- if (identical(position, as.numeric(seq_len(p)))) "Predictor" else "Wavelength"
      return(list(position = position, label = label))
    }
  }
  list(position = as.numeric(seq_len(p)), label = "Predictor")
}
