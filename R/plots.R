#' Plots of LW-PLS local models
#'
#' @description
#' Publication-ready plots of the local models returned by [lwpls_local()]:
#'
#' * `plot_reliability()`: the \eqn{T^2} and \eqn{Q} statistics of the
#'   queries divided by their theoretical or empirical limits, or their
#'   percentile ranks among the training samples of their local models (see
#'   [lwpls_local()]). Queries beyond a limit are extrapolations of their
#'   local model, and their predictions should be checked.
#' * `plot_coefficients()`: the regression coefficients (or VIP) of the local
#'   models along the spectral axis, as one line per query colored by its
#'   prediction, or as a heatmap of the queries ordered by their prediction.
#'   They show which bands drive the predictions and how they change across
#'   the range of the outcome.
#' * `plot_robust_weights()`: the robust weights of the training samples
#'   (robust models), averaged over the local models in which they take part.
#'   Samples with low weights, often wrong reference values, are labeled.
#' * `plot_selection()`: how often each predictor is selected by the sparse
#'   local models (sparse models): for each component (default), or by the
#'   local models as a whole (any component), with the mean spectrum for
#'   reference. On smooth spectra, the components of a sparse model often
#'   select different bands, so that almost all predictors are used by the
#'   model as a whole: the selection per component is then more informative.
#'
#' All plots are ggplot2 objects that can be modified further.
#'
#' @param object An `lwpls_local` object from [lwpls_local()].
#' @param type For `autoplot()`, the plot: `"reliability"`, `"coefficients"`,
#'   `"robust_weights"` or `"selection"`. For `plot_coefficients()`,
#'   `"coefficients"` (default) or `"vip"`.
#' @param labels The number of points to label with their row number: the
#'   queries with the largest statistics beyond the limits, or the training
#'   samples with the smallest robust weights.
#' @param reference For `plot_reliability()`: `"empirical"` or
#'   `"theoretical"` to plot the statistics divided by these limits, or
#'   `"percentile"` to plot their percentile ranks. Defaults to the `limits`
#'   used in [lwpls_local()].
#' @param log A logical: should the axes be on a log scale? Not available for
#'   percentile ranks.
#' @param style `"lines"` (default) or `"heatmap"`.
#' @param outcome The outcome (or class) whose coefficients are shown, or that
#'   is used on the x axis of `plot_robust_weights()`. Defaults to the first.
#' @param reverse A logical: should the spectral axis be reversed, as is usual
#'   for wavenumbers?
#' @param axis_label The label of the spectral axis. Defaults to
#'   `"Wavelength"` or `"Predictor"`, see [lwpls_local()].
#' @param by For `plot_selection()`, `"component"` (default) for the selection
#'   frequency of each predictor in each component, or `"model"` for the
#'   fraction of local models that use each predictor in any component.
#' @param spectrum A logical: should the mean spectrum of the training data,
#'   rescaled to the range of the selection frequencies, be drawn? Only for
#'   `by = "model"`.
#' @param ... Arguments passed to the plot function, for `autoplot()`.
#'   Not used otherwise.
#' @return A ggplot object.
#' @examples
#' fit <- lwpls_fit(mpg ~ ., data = mtcars[-(1:5), ], num_comp = 3, robust = TRUE, sparsity = 0.3)
#' local <- lwpls_local(fit, mtcars[1:5, ])
#'
#' plot_reliability(local)
#' plot_coefficients(local)
#' plot_robust_weights(local)
#' plot_selection(local)
#' @name lwpls_plots
NULL

#' @rdname lwpls_plots
#' @export
plot_reliability <- function(object, reference = NULL, labels = 5, log = FALSE, ...) {
  check_local(object)
  reference <- if (is.null(reference)) object$info$limits else reference
  reference <- rlang::arg_match0(reference, c("empirical", "theoretical", "percentile"))
  check_whole(labels, "labels", min = 0)
  check_bool(log, "log")
  rlang::check_dots_empty()
  if (log && reference == "percentile") {
    cli::cli_abort("{.code log = TRUE} is not available for percentile ranks.")
  }

  rel <- object$reliability
  level <- object$info$level
  if (reference == "percentile") {
    # Equivalent to the empirical limits, on a continuous scale.
    rel$status <- reliability_status(rel, reliability_ratios(rel, "empirical"))
    rel$x <- rel$t2_percentile
    rel$y <- rel$q_percentile
    cutoff <- level
  } else {
    ratios <- reliability_ratios(rel, reference)
    rel$status <- reliability_status(rel, ratios)
    rel$x <- ratios$t2
    rel$y <- ratios$q
    cutoff <- 1
  }

  hidden <- sum(rel$status %in% "few neighbors")
  data <- rel[!is.na(rel$x) & !is.na(rel$y), , drop = FALSE]
  if (log) {
    data <- data[data$x > 0 & data$y > 0, , drop = FALSE]
  }
  flagged <- data[data$status != "inside", , drop = FALSE]
  extreme <- if (reference == "percentile") {
    pmax(rel$t2_ratio, rel$q_ratio)[match(flagged$.row, rel$.row)]
  } else {
    pmax(flagged$x, flagged$y)
  }
  flagged <- flagged[order(extreme, decreasing = TRUE)[seq_len(min(labels, nrow(flagged)))], , drop = FALSE]

  percent <- format(100 * level)
  if (reference == "percentile") {
    x_label <- expression(italic(T)^2 ~ "percentile rank")
    y_label <- expression(italic(Q) ~ "percentile rank")
  } else {
    kind <- if (reference == "empirical") "empirical" else "theoretical"
    x_label <- bquote(italic(T)^2 ~ "/" ~ .(paste0(percent, "% ", kind, " limit")))
    y_label <- bquote(italic(Q) ~ "/" ~ .(paste0(percent, "% ", kind, " limit")))
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(.data$x, .data$y, colour = .data$status, shape = .data$status)
  ) +
    ggplot2::geom_vline(xintercept = cutoff, linetype = 2, colour = "grey50") +
    ggplot2::geom_hline(yintercept = cutoff, linetype = 2, colour = "grey50") +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::geom_text(
      data = flagged,
      ggplot2::aes(label = .data$.row),
      vjust = -0.8,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = status_colours) +
    ggplot2::scale_shape_manual(values = status_shapes) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      colour = NULL,
      shape = NULL,
      caption = if (hidden > 0) {
        paste(hidden, if (hidden == 1) "query with too few neighbors is" else "queries with too few neighbors are", "not shown")
      }
    ) +
    theme_lwpls()
  if (reference == "percentile") {
    plot <- plot +
      ggplot2::scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1))
  } else if (log) {
    plot <- plot + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }
  plot
}

#' @rdname lwpls_plots
#' @export
plot_coefficients <- function(object,
                              type = c("coefficients", "vip"),
                              style = c("lines", "heatmap"),
                              outcome = NULL,
                              reverse = FALSE,
                              axis_label = NULL,
                              ...) {
  check_local(object)
  type <- rlang::arg_match(type)
  style <- rlang::arg_match(style)
  check_bool(reverse, "reverse")
  rlang::check_dots_empty()
  outcome <- check_local_outcome(object, outcome)

  coefs <- object$coefficients
  coefs <- coefs[coefs$outcome == outcome & !is.na(coefs$coefficient), , drop = FALSE]
  coefs$value <- if (type == "vip") coefs$vip else coefs$coefficient
  value_label <- if (type == "vip") "VIP" else "Coefficient"

  # Prediction of each query, to color the lines and order the heatmap.
  rel <- object$reliability
  regression <- object$info$mode == "regression"
  pred_col <- if (!regression) {
    ".pred_class"
  } else if (length(object$info$outcome_names) == 1) {
    ".pred"
  } else {
    paste0(".pred_", outcome)
  }
  coefs$prediction <- rel[[pred_col]][match(coefs$.row, rel$.row)]

  if (style == "lines") {
    plot <- ggplot2::ggplot(
      coefs,
      ggplot2::aes(.data$position, .data$value, group = .data$.row, colour = .data$prediction)
    )
    if (type == "coefficients") {
      plot <- plot + ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3)
    }
    plot <- plot +
      ggplot2::geom_line(linewidth = 0.4, alpha = 0.8) +
      ggplot2::labs(colour = if (regression) "Prediction" else "Predicted class")
    plot <- plot + if (regression) {
      ggplot2::scale_colour_viridis_c(option = "mako", end = 0.9)
    } else {
      ggplot2::scale_colour_viridis_d(option = "mako", end = 0.9)
    }
    y_label <- value_label
  } else {
    queries <- unique(coefs[c(".row", "prediction")])
    queries <- queries[order(queries$prediction, queries$.row), , drop = FALSE]
    coefs$query <- match(coefs$.row, queries$.row)
    plot <- ggplot2::ggplot(coefs, ggplot2::aes(.data$position, .data$query, fill = .data$value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(fill = value_label)
    plot <- plot + if (type == "vip") {
      ggplot2::scale_fill_viridis_c(option = "mako", direction = -1)
    } else {
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B")
    }
    y_label <- "Queries, ordered by prediction"
  }

  plot <- plot +
    ggplot2::labs(
      x = if (is.null(axis_label)) object$info$axis_label else axis_label,
      y = y_label,
      subtitle = if (!regression || length(object$info$outcome_names) > 1) paste("Outcome:", outcome)
    ) +
    theme_lwpls()
  if (style == "heatmap") {
    plot <- plot + ggplot2::scale_x_continuous(expand = c(0, 0), transform = if (reverse) "reverse" else "identity")
  } else if (reverse) {
    plot <- plot + ggplot2::scale_x_reverse()
  }
  plot
}

#' @rdname lwpls_plots
#' @export
plot_robust_weights <- function(object, labels = 5, outcome = NULL, ...) {
  check_local(object)
  check_whole(labels, "labels", min = 0)
  rlang::check_dots_empty()
  if (!object$info$robust) {
    cli::cli_abort(c(
      "Robust weights are only available for robust models.",
      "i" = "Fit the model with {.code robust = TRUE}."
    ))
  }

  training <- object$training[!is.na(object$training$robust_weight), , drop = FALSE]
  regression <- object$info$mode == "regression"
  if (regression) {
    outcome <- check_local_outcome(object, outcome)
    training$x <- training[[outcome]]
    x_label <- if (outcome == ".outcome") "Reference value" else paste0("Reference value (", outcome, ")")
  } else {
    training$x <- training$.class
    x_label <- "Class"
  }
  lowest <- training[order(training$robust_weight)[seq_len(min(labels, nrow(training)))], , drop = FALSE]

  point <- if (regression) {
    ggplot2::geom_point(ggplot2::aes(size = .data$influence), colour = "#08519C", alpha = 0.7)
  } else {
    ggplot2::geom_jitter(
      ggplot2::aes(size = .data$influence),
      colour = "#08519C",
      alpha = 0.7,
      width = 0.15,
      height = 0
    )
  }
  ggplot2::ggplot(training, ggplot2::aes(.data$x, .data$robust_weight)) +
    point +
    ggplot2::geom_text(
      data = lowest,
      ggplot2::aes(label = .data$.row),
      vjust = -0.9,
      size = 3,
      colour = "#B2182B",
      check_overlap = TRUE
    ) +
    ggplot2::scale_size_area(max_size = 3.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(x = x_label, y = "Robust weight", size = "Influence") +
    theme_lwpls()
}

#' @rdname lwpls_plots
#' @export
plot_selection <- function(object,
                           by = c("component", "model"),
                           spectrum = TRUE,
                           reverse = FALSE,
                           axis_label = NULL,
                           ...) {
  check_local(object)
  by <- rlang::arg_match(by)
  check_bool(spectrum, "spectrum")
  check_bool(reverse, "reverse")
  rlang::check_dots_empty()
  if (object$info$sparsity <= 0) {
    cli::cli_abort(c(
      "Predictor selection is only available for sparse models.",
      "i" = "Fit the model with {.code sparsity > 0}."
    ))
  }
  x_label <- if (is.null(axis_label)) object$info$axis_label else axis_label

  if (by == "component") {
    data <- object$selection
    plot <- ggplot2::ggplot(data, ggplot2::aes(.data$position, .data$component, fill = .data$frequency)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(
        option = "mako",
        direction = -1,
        limits = c(0, 1),
        labels = scales::label_percent()
      ) +
      ggplot2::scale_y_reverse(breaks = unique(data$component), expand = c(0, 0)) +
      ggplot2::labs(x = x_label, y = "Component", fill = "Selected") +
      theme_lwpls() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    return(plot + ggplot2::scale_x_continuous(
      expand = c(0, 0),
      transform = if (reverse) "reverse" else "identity"
    ))
  }

  coefs <- object$coefficients
  coefs <- coefs[coefs$outcome == object$info$outcome_names[1] & !is.na(coefs$selected), , drop = FALSE]
  frequency <- tapply(coefs$selected, coefs$position, mean)
  data <- tibble::tibble(position = as.numeric(names(frequency)), frequency = as.vector(frequency))

  bar_width <- if (nrow(data) > 1) min(diff(sort(data$position))) else 1
  plot <- ggplot2::ggplot(data, ggplot2::aes(.data$position, .data$frequency)) +
    ggplot2::geom_col(fill = "#08519C", width = bar_width, alpha = 0.85)
  if (spectrum) {
    center <- object$info$center
    rng <- range(center)
    scaled <- if (diff(rng) > 0) (center - rng[1]) / diff(rng) else rep(0.5, length(center))
    spec <- tibble::tibble(position = object$info$position, value = scaled)
    plot <- plot + ggplot2::geom_line(
      data = spec,
      ggplot2::aes(.data$position, .data$value),
      colour = "grey35",
      linewidth = 0.5
    )
  }
  plot <- plot +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    ggplot2::labs(
      x = x_label,
      y = "Local models using the predictor",
      subtitle = if (spectrum) "Grey line: mean spectrum (rescaled)"
    ) +
    theme_lwpls()
  if (reverse) {
    plot <- plot + ggplot2::scale_x_reverse()
  }
  plot
}

#' @rdname lwpls_plots
#' @export
autoplot.lwpls_local <- function(object,
                                 type = c("reliability", "coefficients", "robust_weights", "selection"),
                                 ...) {
  type <- rlang::arg_match(type)
  switch(
    type,
    reliability = plot_reliability(object, ...),
    coefficients = plot_coefficients(object, ...),
    robust_weights = plot_robust_weights(object, ...),
    selection = plot_selection(object, ...)
  )
}

# ------------------------------------------------------------------------------

status_colours <- c(
  "inside" = "#4D4D4D",
  "high T2" = "#E69F00",
  "high Q" = "#0072B2",
  "high T2 and Q" = "#D55E00",
  "few neighbors" = "#CC79A7"
)

status_shapes <- c(
  "inside" = 16,
  "high T2" = 17,
  "high Q" = 15,
  "high T2 and Q" = 18,
  "few neighbors" = 4
)

theme_lwpls <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      plot.title.position = "plot",
      legend.key = ggplot2::element_blank()
    )
}

check_local <- function(object, call = rlang::caller_env()) {
  if (!inherits(object, "lwpls_local")) {
    cli::cli_abort(
      "{.arg object} must be the result of {.fn lwpls_local}, not {.obj_type_friendly {object}}.",
      call = call
    )
  }
  invisible(object)
}

check_local_outcome <- function(object, outcome, call = rlang::caller_env()) {
  values <- object$info$outcome_names
  if (is.null(outcome)) {
    return(values[1])
  }
  if (!rlang::is_string(outcome) || !outcome %in% values) {
    cli::cli_abort(
      "{.arg outcome} must be one of {.or {.val {values}}}, not {describe(outcome)}.",
      call = call
    )
  }
  outcome
}
