#' Localization parameter for LW-PLS
#'
#' The localization parameter controls how fast the similarity weights of
#' [lwpls()] decrease with the distance to the query sample. Small values
#' give very local models; large values give all training samples similar
#' weights, so that the model tends to a global PLS model.
#'
#' The default range, \eqn{2^{-9}}{2^-9} to \eqn{2^{5}}{2^5} on a log-2 scale, follows
#' the candidate values used by Kaneko for the LW-PLS hyperparameter search.
#'
#' @param range A two-element vector holding the _defaults_ for the smallest
#'   and largest possible values, respectively. If a transformation is
#'   specified, these values should be in the _transformed units_.
#' @param trans A `trans` object from the `scales` package, such as
#'   [scales::transform_log10()] or [scales::transform_reciprocal()]. If not
#'   provided, the default is used which matches the units used in `range`.
#'   If no transformation, `NULL`.
#' @return A `quant_param` object (see [dials::new_quant_param()]).
#' @examples
#' localization()
#' dials::value_seq(localization(), 5)
#' dials::grid_regular(dials::num_comp(c(1, 5)), localization(), levels = 3)
#' @export
localization <- function(range = c(-9, 5), trans = scales::transform_log2()) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    trans = trans,
    label = c(localization = "Localization Parameter"),
    finalize = NULL
  )
}
