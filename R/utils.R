# Argument checkers ------------------------------------------------------------

check_whole <- function(x,
                        arg,
                        min = -Inf,
                        allow_null = FALSE,
                        call = rlang::caller_env()) {
  if (is.null(x) && allow_null) {
    return(invisible(NULL))
  }
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok) {
    cli::cli_abort(
      "{.arg {arg}} must be a single whole number >= {min}, not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_positive <- function(x, arg, allow_null = FALSE, call = rlang::caller_env()) {
  if (is.null(x) && allow_null) {
    return(invisible(NULL))
  }
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) && x > 0
  if (!ok) {
    cli::cli_abort(
      "{.arg {arg}} must be a single positive number, not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_similarity <- function(x, allow_null = FALSE, call = rlang::caller_env()) {
  if (is.null(x) && allow_null) {
    return(invisible(NULL))
  }
  if (!rlang::is_string(x) || !x %in% values_similarity) {
    cli::cli_abort(
      "{.arg similarity} must be one of {.or {.val {values_similarity}}}, not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_sparsity <- function(x, allow_null = FALSE, call = rlang::caller_env()) {
  if (is.null(x) && allow_null) {
    return(invisible(NULL))
  }
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 0 && x < 1
  if (!ok) {
    cli::cli_abort(
      "{.arg sparsity} must be a single number in [0, 1), not {describe(x)}.",
      call = call
    )
  }
  invisible(x)
}

check_weight_function <- function(x, call = rlang::caller_env()) {
  values <- c("fair", "hampel")
  if (!rlang::is_string(x) || !x %in% values) {
    cli::cli_abort(
      "{.arg weight_function} must be one of {.or {.val {values}}}, not {describe(x)}.",
      call = call
    )
  }
  x
}

# Returns the tuning constant(s), with the defaults for NULL.
check_robust_constant <- function(x, weight_function, call = rlang::caller_env()) {
  if (weight_function == "fair") {
    if (is.null(x)) {
      return(4)
    }
    check_positive(x, "robust_constant", call = call)
    return(as.numeric(x))
  }
  if (is.null(x)) {
    return(c(0.95, 0.975, 0.999))
  }
  ok <- is.numeric(x) && length(x) == 3 && !anyNA(x) && all(x > 0 & x < 1) &&
    all(diff(x) > 0)
  if (!ok) {
    cli::cli_abort(
      c(
        "{.arg robust_constant} must be three increasing probabilities for the
         Hampel function, not {describe(x)}.",
        "i" = "For example, {.code c(0.95, 0.975, 0.999)}."
      ),
      call = call
    )
  }
  as.numeric(x)
}

check_bool <- function(x, arg, call = rlang::caller_env()) {
  if (!rlang::is_bool(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be {.code TRUE} or {.code FALSE}, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  invisible(x)
}

# Short description of an invalid argument value for error messages.
describe <- function(x) {
  if (is.numeric(x) && length(x) == 1) {
    format(x)
  } else if (rlang::is_string(x)) {
    encodeString(x, quote = '"')
  } else {
    cli::format_inline("{.obj_type_friendly {x}}")
  }
}
