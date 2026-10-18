#' Locally-weighted partial least squares
#'
#' @description
#' `lwpls()` defines a locally-weighted partial least squares (LW-PLS) model
#' for use with the [parsnip][parsnip::parsnip] package. LW-PLS is a
#' *just-in-time* (lazy) learner: instead of fitting one global model, a
#' dedicated weighted PLS model is built around every new sample, where each
#' training sample is weighted by its similarity to that sample. This lets a
#' linear method capture nonlinear and time-varying relationships, which makes
#' it popular for soft sensors and spectroscopic calibration.
#'
#' The model can be used for regression (one or more numeric outcomes) and
#' classification (discriminant analysis on class indicators). The local
#' models can also be sparse, so that each one selects its own predictors,
#' and robust to outliers in the training data.
#'
#' There is a single engine, `"lwpls"`, implemented in this package (C++ via
#' RcppArmadillo).
#'
#' @param mode A single character string for the prediction outcome mode.
#'   Possible values are `"unknown"`, `"regression"`, or `"classification"`.
#' @param num_comp The number of PLS components (latent variables) in each
#'   local model (engine default: 2).
#' @param localization A positive number controlling how local the models
#'   are (engine default: 1). Small values give large weights to the closest
#'   training samples only; large values give all samples similar weights, so
#'   that LW-PLS tends to a global PLS model. See Details.
#' @param neighbors The number of nearest training samples used for each
#'   local model. The default (`NULL`) uses all training samples, as in the
#'   original LW-PLS algorithm.
#' @param similarity The similarity index: `"euclidean"` (engine default) or
#'   `"covariance"` for covariance-based LW-PLS (CbLW-PLS). See Details and
#'   [similarity()].
#' @param sparsity A number in `[0, 1)`: the sparsity threshold \eqn{\eta} of
#'   the local models (engine default: 0, no sparsity). See Details and
#'   [sparsity()].
#' @param engine A single character string specifying the computational
#'   engine. Only `"lwpls"` is available.
#'
#' @details
#' For a query sample \eqn{x_q}, the similarity weight of the training sample
#' \eqn{x_i} is
#'
#' \deqn{\omega_i = \exp\left(-\frac{d_i}{\sigma_d \, \varphi}\right),}{w_i = exp(-d_i / (sigma_d * phi)),}
#'
#' where \eqn{d_i = \lVert x_i - x_q \rVert}{d_i = ||x_i - x_q||} is the Euclidean distance (on
#' standardized predictors by default), \eqn{\sigma_d} is the standard
#' deviation of the distances and \eqn{\varphi}{phi} is the `localization`
#' parameter (\eqn{\lambda} in Kaneko's implementation). A weighted PLS model
#' is then fitted with weighted centering, weighted covariances and the
#' query is projected onto its local latent space to obtain the prediction.
#'
#' With `similarity = "covariance"`, the distances are computed after
#' projecting the samples on the covariance direction
#' \eqn{\Gamma = X^\top Y / \lVert X^\top Y \rVert}{Gamma = X'Y / ||X'Y||} of the training data,
#' \eqn{d_i = \lVert \Gamma^\top (x_i - x_q) \rVert}{d_i = ||Gamma'(x_i - x_q)||}. This covariance-based
#' LW-PLS (CbLW-PLS; Hazama and Kano, 2015) accounts for the relationships
#' among the predictors and between the predictors and the outcome(s). For a
#' single outcome, \eqn{\Gamma} is the first PLS weight vector; for several
#' outcomes (or classes), the distance combines the covariance directions of
#' all of them. Hazama and Kano write the weights as
#' \eqn{\exp(-\phi d_i / \sigma_d)}{exp(-phi d_i / sigma_d)}, so their \eqn{\phi} is
#' `1 / localization`.
#'
#' # Sparse local models
#'
#' With `sparsity` \eqn{\eta > 0}, the local models are fitted by sparse
#' NIPALS (SNIPLS; Hoffmann et al., 2015): the weight vector of each component
#' is soft-thresholded at \eqn{\eta \max_j |w_j|}{eta * max|w|}, so that the
#' predictors with small weights are dropped, and the X loadings are zero
#' outside the predictors selected so far. Each local model therefore selects
#' its own predictors, such as the wavelengths relevant around the query. For
#' several outcomes or classes, the weight vector is the dominant direction of
#' the covariance with all of them before thresholding.
#'
#' # Robust local models
#'
#' With the engine argument `robust = TRUE`, each local model is fitted by
#' partial robust M-regression (PRM; Serneels et al., 2005), which makes it
#' robust to outlying reference values (vertical outliers) and to outlying
#' samples in the predictor space (leverage points) among the neighbors of the
#' query. Each training sample gets the weight
#' \deqn{\omega_i \, w_i^r \, w_i^t,}{w_i * wr_i * wt_i,}
#' the product of its similarity weight \eqn{\omega_i}{w_i}, a weight
#' \eqn{w_i^r}{wr_i} that decreases with its residual in the local model, and a
#' weight \eqn{w_i^t}{wt_i} that decreases with its distance to the center of
#' the scores, and the local model is refitted with the new weights until the
#' norm of its regression coefficients changes by less than 1%. The weights
#' are given by the Fair function
#' \eqn{f(z) = 1 / (1 + |z / c|)^2}{f(z) = 1 / (1 + |z / c|)^2} (default,
#' \eqn{c = 4}) of the residuals and distances divided by their median, or by
#' the Hampel function, which gives zero weight to the most outlying samples,
#' with cutoffs at quantiles of the \eqn{\chi}{chi} distribution (Hoffmann et
#' al., 2015). As in PRM, the local center of the predictors is their
#' L1-median, the intercept is the median of the residuals, and all these
#' medians are weighted by the similarity weights, so the model stays local.
#' For several outcomes, the residual weights use the norm of the scaled
#' residual vectors; for classification, the indicators are centered by their
#' (weighted) means so that the predicted indicators still sum to one.
#'
#' The robust models also standardize the data by the median and the median
#' absolute deviation (MAD), and use the MAD of the distances instead of their
#' standard deviation in the similarity weights. With
#' `similarity = "covariance"`, \eqn{\Gamma}{Gamma} is computed with the
#' initial (Fair) weights of PRM.
#'
#' The weights depend on the number of components, so a robust model is
#' refitted for each value of `num_comp`: tuning `num_comp` costs about as
#' many robust fits as values in the grid.
#'
#' # Neighbors
#'
#' When `neighbors` is set, only the `neighbors` closest training samples
#' receive a non-zero weight and \eqn{\sigma_d} is computed over their
#' distances. This is the "KNN-LW" strategy, which Lesnoff et al. (2020)
#' compare with the original LW-PLS that weights all training samples.
#'
#' For classification, the outcome is converted to class indicator columns
#' and modeled with a multi-response LW-PLS (LW-PLS-DA; Bevilacqua and Marini,
#' 2014). The predicted class
#' is the one with the largest predicted indicator. Class probabilities are
#' obtained by truncating the predicted indicators to `[0, 1]` and
#' renormalizing them to sum to one.
#'
#' Because LW-PLS is a lazy learner, fitting only stores the (standardized)
#' training data; the computational work happens at prediction time. A
#' prediction for `num_comp` components also yields the predictions for all
#' smaller numbers of components, so [tune::tune_grid()] evaluates every value
#' of `num_comp` from a single model (the "submodel trick").
#'
#' # Engine details
#'
#' The `"lwpls"` engine calls [lwpls_fit()]. It accepts the engine arguments:
#'
#' * `scale` (logical, default `TRUE`): should the predictors (and numeric
#'   outcomes) be standardized to unit variance? Set it to `FALSE` for spectra
#'   and other data measured on a common scale.
#' * `robust` (logical, default `FALSE`): robust local models (see above).
#' * `weight_function` (`"fair"`, the default, or `"hampel"`),
#'   `robust_constant` and `max_iter` (default 30): the settings of the robust
#'   models, see [lwpls_fit()].
#'
#' Predictors with factor columns are converted to dummy variables when the
#' model is fitted with a formula.
#'
#' Tuning parameters:
#'
#' * `num_comp`: [dials::num_comp()]
#' * `localization`: [localization()]
#' * `neighbors`: [dials::neighbors()] (default range 10--200)
#' * `similarity`: [similarity()]
#' * `sparsity`: [sparsity()]
#'
#' @references
#' Kim, S., Kano, M., Nakagawa, H. and Hasebe, S. (2011). Estimation of active
#' pharmaceutical ingredients content using locally weighted partial least
#' squares and statistical wavelength selection. *International Journal of
#' Pharmaceutics*, 421(2), 269--274. \doi{10.1016/j.ijpharm.2011.10.007}
#'
#' Bevilacqua, M. and Marini, F. (2014). Local classification: Locally
#' weighted-partial least squares-discriminant analysis (LW-PLS-DA).
#' *Analytica Chimica Acta*, 838, 20--30. \doi{10.1016/j.aca.2014.05.057}
#'
#' Hazama, K. and Kano, M. (2015). Covariance-based locally weighted partial
#' least squares for high-performance adaptive modeling. *Chemometrics and
#' Intelligent Laboratory Systems*, 146, 55--62.
#' \doi{10.1016/j.chemolab.2015.05.007}
#'
#' Serneels, S., Croux, C., Filzmoser, P. and Van Espen, P. J. (2005). Partial
#' robust M-regression. *Chemometrics and Intelligent Laboratory Systems*,
#' 79(1--2), 55--64. \doi{10.1016/j.chemolab.2005.04.007}
#'
#' Hoffmann, I., Serneels, S., Filzmoser, P. and Croux, C. (2015). Sparse
#' partial robust M regression. *Chemometrics and Intelligent Laboratory
#' Systems*, 149, 50--59. \doi{10.1016/j.chemolab.2015.09.019}
#'
#' Lesnoff, M., Metz, M. and Roger, J.-M. (2020). Comparison of locally
#' weighted PLS strategies for regression and discrimination on agronomic NIR
#' data. *Journal of Chemometrics*, 34(5), e3209. \doi{10.1002/cem.3209}
#'
#' @return A model specification object with classes `lwpls` and
#'   `model_spec`.
#' @seealso [lwpls_fit()] for the underlying fitting function,
#'   [multi_predict._lwpls_fit()], [localization()], [similarity()],
#'   [sparsity()].
#' @examples
#' lwpls(num_comp = 3, localization = 0.5) |>
#'   parsnip::set_mode("regression") |>
#'   parsnip::translate()
#'
#' lwpls_mod <-
#'   lwpls(num_comp = 3, localization = 0.5) |>
#'   parsnip::set_mode("regression") |>
#'   parsnip::fit(mpg ~ ., data = mtcars[-(1:5), ])
#'
#' lwpls_mod
#' predict(lwpls_mod, new_data = mtcars[1:5, ])
#'
#' # Robust and sparse local models
#' lwpls(num_comp = 3, sparsity = 0.3) |>
#'   parsnip::set_mode("regression") |>
#'   parsnip::set_engine("lwpls", robust = TRUE) |>
#'   parsnip::fit(mpg ~ ., data = mtcars[-(1:5), ]) |>
#'   predict(new_data = mtcars[1:5, ])
#'
#' # Classification
#' lwpls(num_comp = 2) |>
#'   parsnip::set_mode("classification") |>
#'   parsnip::fit(Species ~ ., data = iris) |>
#'   predict(new_data = iris[c(1, 51, 101), ], type = "prob")
#' @export
lwpls <- function(mode = "unknown",
                  num_comp = NULL,
                  localization = NULL,
                  neighbors = NULL,
                  similarity = NULL,
                  sparsity = NULL,
                  engine = "lwpls") {
  args <- list(
    num_comp = rlang::enquo(num_comp),
    localization = rlang::enquo(localization),
    neighbors = rlang::enquo(neighbors),
    similarity = rlang::enquo(similarity),
    sparsity = rlang::enquo(sparsity)
  )

  parsnip::new_model_spec(
    "lwpls",
    args = args,
    eng_args = NULL,
    mode = mode,
    user_specified_mode = !missing(mode),
    method = NULL,
    engine = engine,
    user_specified_engine = !missing(engine)
  )
}

#' @export
print.lwpls <- function(x, ...) {
  parsnip::print_model_spec(x, desc = "Locally-Weighted PLS", ...)
}

#' Update a locally-weighted PLS specification
#'
#' @param object A [lwpls()] model specification.
#' @param parameters A 1-row tibble or named list with *main* parameters to
#'   update. Use **either** `parameters` **or** the main arguments directly
#'   when updating. If the main arguments are used, these will supersede the
#'   values in `parameters`. Also, using engine arguments in this object will
#'   result in an error.
#' @inheritParams lwpls
#' @param fresh A logical for whether the arguments should be modified
#'   in-place or replaced wholesale.
#' @param ... Not used for `update()`.
#' @return An updated model specification.
#' @examples
#' spec <- lwpls(num_comp = 3)
#' spec
#' update(spec, localization = 0.25)
#' update(spec, localization = 0.25, fresh = TRUE)
#' @export
update.lwpls <- function(object,
                         parameters = NULL,
                         num_comp = NULL,
                         localization = NULL,
                         neighbors = NULL,
                         similarity = NULL,
                         sparsity = NULL,
                         fresh = FALSE,
                         ...) {
  args <- list(
    num_comp = rlang::enquo(num_comp),
    localization = rlang::enquo(localization),
    neighbors = rlang::enquo(neighbors),
    similarity = rlang::enquo(similarity),
    sparsity = rlang::enquo(sparsity)
  )

  parsnip::update_spec(
    object = object,
    parameters = parameters,
    args_enquo_list = args,
    fresh = fresh,
    cls = "lwpls",
    ...
  )
}

#' @exportS3Method parsnip::check_args
check_args.lwpls <- function(object, call = rlang::caller_env()) {
  args <- lapply(object$args, rlang::eval_tidy)
  check_whole(args$num_comp, "num_comp", min = 1, allow_null = TRUE, call = call)
  check_positive(args$localization, "localization", allow_null = TRUE, call = call)
  check_whole(args$neighbors, "neighbors", min = 1, allow_null = TRUE, call = call)
  check_similarity(args$similarity, allow_null = TRUE, call = call)
  check_sparsity(args$sparsity, allow_null = TRUE, call = call)
  invisible(object)
}

#' Determine the minimal set of model fits for tuning
#'
#' `min_grid()` is used by the tune package to fit as few models as possible.
#' LW-PLS predictions for `num_comp` components include those for all
#' smaller numbers of components, so only the largest `num_comp` value of
#' each combination of the other parameters is fitted.
#'
#' @param x A [lwpls()] model specification.
#' @param grid A tibble with tuning parameter combinations.
#' @param ... Not currently used.
#' @return A tibble with the minimal tuning grid and a list column with the
#'   submodel values to predict.
#' @keywords internal
#' @export
min_grid.lwpls <- function(x, grid, ...) {
  rlang::check_installed("tune")
  tune::fit_max_value(x, grid, ...)
}
