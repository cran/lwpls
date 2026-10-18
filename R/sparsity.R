#' Sparsity of the LW-PLS local models
#'
#' The sparsity threshold \eqn{\eta} of [lwpls()]: the weight vector of each
#' component of a local model is soft-thresholded at
#' \eqn{\eta \max_j |w_j|}{eta * max|w|} (sparse NIPALS, Hoffmann et al.,
#' 2015). Zero gives ordinary local PLS models; larger values drop more
#' predictors.
#'
#' @param range A two-element vector holding the _defaults_ for the smallest
#'   and largest possible values, respectively.
#' @param trans A `trans` object from the `scales` package, or `NULL` (the
#'   default) for no transformation.
#' @return A `quant_param` object (see [dials::new_quant_param()]).
#' @references
#' Hoffmann, I., Serneels, S., Filzmoser, P. and Croux, C. (2015). Sparse
#' partial robust M regression. *Chemometrics and Intelligent Laboratory
#' Systems*, 149, 50--59. \doi{10.1016/j.chemolab.2015.09.019}
#' @examples
#' sparsity()
#' dials::value_seq(sparsity(), 4)
#' @export
sparsity <- function(range = c(0, 0.9), trans = NULL) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    trans = trans,
    label = c(sparsity = "Sparsity Threshold"),
    finalize = NULL
  )
}
