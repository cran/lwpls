#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib lwpls, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom generics min_grid
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#' @importFrom parsnip multi_predict
#' @importFrom stats predict update
## usethis namespace: end
NULL

# The model is registered in parsnip's model environment when the package is
# loaded. The check avoids re-registering it when the namespace is reloaded
# (e.g. by `devtools::load_all()`).
# nocov start
.onLoad <- function(libname, pkgname) {
  current <- parsnip::get_model_env()
  if (!any(current$models == "lwpls")) {
    make_lwpls()
  }
}
# nocov end
