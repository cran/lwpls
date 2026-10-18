#' Similarity index for LW-PLS
#'
#' The similarity index determines how the distance between a query sample
#' and the training samples is measured before it is turned into weights by
#' [lwpls()]:
#'
#' * `"euclidean"`: the Euclidean distance between the (standardized)
#'   predictors, as in the original LW-PLS.
#' * `"covariance"`: the distance after projecting the samples on the
#'   covariance direction \eqn{\Gamma = X^\top Y / \lVert X^\top Y \rVert}{Gamma = X'Y / ||X'Y||}
#'   (covariance-based LW-PLS, CbLW-PLS; Hazama and Kano, 2015). It accounts
#'   for the relationships among the predictors and between the predictors and
#'   the outcome(s).
#'
#' See `vignette("similarity")` for how to choose between them.
#'
#' @param values A character vector of possible similarity indexes.
#' @return A `qual_param` object (see [dials::new_qual_param()]).
#' @references
#' Hazama, K. and Kano, M. (2015). Covariance-based locally weighted partial
#' least squares for high-performance adaptive modeling. *Chemometrics and
#' Intelligent Laboratory Systems*, 146, 55--62.
#' \doi{10.1016/j.chemolab.2015.05.007}
#' @examples
#' similarity()
#' values_similarity
#' @export
similarity <- function(values = values_similarity) {
  dials::new_qual_param(
    type = "character",
    values = values,
    label = c(similarity = "Similarity Index"),
    finalize = NULL
  )
}

#' @rdname similarity
#' @export
values_similarity <- c("euclidean", "covariance")
