# Registration of the `lwpls` model and its `"lwpls"` engine in parsnip's
# model environment. Called from `.onLoad()`, i.e. before covr can trace it.

# nocov start
make_lwpls <- function() {
  parsnip::set_new_model("lwpls")

  for (mode in c("regression", "classification")) {
    parsnip::set_model_mode(model = "lwpls", mode = mode)
    parsnip::set_model_engine(model = "lwpls", mode = mode, eng = "lwpls")
    parsnip::set_dependency(model = "lwpls", eng = "lwpls", pkg = "lwpls", mode = mode)

    parsnip::set_fit(
      model = "lwpls",
      eng = "lwpls",
      mode = mode,
      value = list(
        interface = "matrix",
        protect = c("x", "y"),
        func = c(pkg = "lwpls", fun = "lwpls_fit"),
        defaults = list()
      )
    )

    parsnip::set_encoding(
      model = "lwpls",
      eng = "lwpls",
      mode = mode,
      options = list(
        predictor_indicators = "traditional",
        compute_intercept = FALSE,
        remove_intercept = TRUE,
        allow_sparse_x = FALSE
      )
    )

    parsnip::set_pred(
      model = "lwpls",
      eng = "lwpls",
      mode = mode,
      type = "raw",
      value = lwpls_pred_value("raw")
    )
  }

  parsnip::set_model_arg(
    model = "lwpls",
    eng = "lwpls",
    parsnip = "num_comp",
    original = "num_comp",
    func = list(pkg = "dials", fun = "num_comp"),
    has_submodel = TRUE
  )
  parsnip::set_model_arg(
    model = "lwpls",
    eng = "lwpls",
    parsnip = "localization",
    original = "localization",
    func = list(pkg = "lwpls", fun = "localization"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "lwpls",
    eng = "lwpls",
    parsnip = "neighbors",
    original = "neighbors",
    func = list(pkg = "dials", fun = "neighbors", range = c(10L, 200L)),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "lwpls",
    eng = "lwpls",
    parsnip = "similarity",
    original = "similarity",
    func = list(pkg = "lwpls", fun = "similarity"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "lwpls",
    eng = "lwpls",
    parsnip = "sparsity",
    original = "sparsity",
    func = list(pkg = "lwpls", fun = "sparsity"),
    has_submodel = FALSE
  )

  parsnip::set_pred(
    model = "lwpls",
    eng = "lwpls",
    mode = "regression",
    type = "numeric",
    value = lwpls_pred_value("numeric")
  )
  parsnip::set_pred(
    model = "lwpls",
    eng = "lwpls",
    mode = "classification",
    type = "class",
    value = lwpls_pred_value("class")
  )
  parsnip::set_pred(
    model = "lwpls",
    eng = "lwpls",
    mode = "classification",
    type = "prob",
    value = lwpls_pred_value("prob")
  )
}

# The engine's predict method already returns parsnip-formatted tibbles.
lwpls_pred_value <- function(type) {
  list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict"),
    args = list(
      object = quote(object$fit),
      new_data = quote(new_data),
      type = type
    )
  )
}
# nocov end
