dat <- sim_reg(n = 80, seed = 3)
x_names <- paste0("x", 1:6)
# Vertical outliers and leverage points in the training data
dat$y[1:6] <- dat$y[1:6] + 8
dat[7:9, x_names] <- dat[7:9, x_names] + 4
train <- dat[1:65, ]
test <- dat[66:80, ]

# Kernels ---------------------------------------------------------------------

test_that("the general kernel reproduces the fast kernel without sparsity", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 4, localization = 0.5)
  new_x <- standardize(as.matrix(test[x_names]), fit$x_center, fit$x_scale)
  fast <- lwpls_predict_cpp(
    fit$x, fit$y, new_x, fit$x, new_x,
    num_comp = 4, localization = 0.5, neighbors = fit$neighbors, tol = lwpls_tol
  )
  general <- lwpls_general_cpp(
    fit$x, fit$y, new_x, fit$x, new_x,
    comps = 1:4, localization = 0.5, neighbors = fit$neighbors, sparsity = 0,
    robust = 0L, fair_c = 4, hampel_probs = c(0.95, 0.975, 0.999), max_iter = 1L,
    classification = FALSE, tol = lwpls_tol
  )
  expect_equal(general, fast, tolerance = 1e-10)
})

# Sparse local models ---------------------------------------------------------

test_that("sparse models match the SNIPLS reference", {
  for (eta in c(0.2, 0.5)) {
    fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, localization = 0.5, sparsity = eta)
    expect_equal(predict(fit, test[x_names])$.pred, ref_predict(fit, test[x_names], 3), tolerance = 1e-10)
  }

  fit2 <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 3, sparsity = 0.3)
  expect_equal(
    unname(as.matrix(predict(fit2, test[x_names]))),
    unname(ref_predict(fit2, test[x_names], 3)),
    tolerance = 1e-10
  )
})

test_that("sparse models give all numbers of components from one fit", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 4, sparsity = 0.3)
  raw <- predict(fit, test[x_names], type = "raw")
  for (a in 1:4) {
    expect_equal(unname(raw[, 1, a]), ref_predict(fit, test[x_names], a), tolerance = 1e-10)
  }
})

test_that("sparse models ignore the predictors they drop", {
  wide <- withr::with_seed(5, matrix(rnorm(60 * 30), 60, 30))
  colnames(wide) <- paste0("v", 1:30)
  y <- wide[, 1] - wide[, 2] + withr::with_seed(6, rnorm(60, sd = 0.1))
  # Global models, so that the similarity weights do not depend on the query
  dense <- lwpls_fit(wide[1:50, ], y[1:50], num_comp = 1, localization = 1e12)
  sparse <- lwpls_fit(wide[1:50, ], y[1:50], num_comp = 1, localization = 1e12, sparsity = 0.5)

  # Predictors dropped by the threshold in the (single) global component
  w <- abs(drop(crossprod(sparse$x, sparse$y)))
  dropped <- which(w < 0.5 * max(w))
  expect_gt(length(dropped), 20)

  query <- wide[51:55, ]
  moved <- query
  moved[, dropped] <- withr::with_seed(7, matrix(rnorm(5 * length(dropped)), 5))
  expect_equal(predict(sparse, moved), predict(sparse, query), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(predict(dense, moved), predict(dense, query))))
})

test_that("sparse classification", {
  fit <- lwpls_fit(Species ~ ., data = iris, num_comp = 2, sparsity = 0.4)
  prob <- predict(fit, iris[c(1, 51, 101), ], type = "prob")
  expect_equal(rowSums(prob), rep(1, 3))
  expect_gt(mean(predict(fit, iris)$.pred_class == iris$Species), 0.9)
})

# Robust local models ---------------------------------------------------------

test_that("robust models match the PRM reference", {
  for (fun in c("fair", "hampel")) {
    for (a in c(1, 3)) {
      fit <- lwpls_fit(
        train[x_names], train$y,
        num_comp = a, localization = 1, robust = TRUE, weight_function = fun
      )
      # Weighted medians can jump between neighboring values when the
      # cumulative weight is within rounding of one half.
      expect_equal(predict(fit, test[x_names])$.pred, ref_predict(fit, test[x_names], a), tolerance = 1e-6)
    }
  }
})

test_that("robust sparse models (SPRM) match the reference", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, robust = TRUE, sparsity = 0.3)
  expect_equal(predict(fit, test[x_names])$.pred, ref_predict(fit, test[x_names], 3), tolerance = 1e-6)
})

test_that("robust models with several outcomes and classes", {
  fit <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 2, robust = TRUE)
  pred <- predict(fit, test[x_names])
  expect_named(pred, c(".pred_y", ".pred_y2"))
  expect_equal(unname(as.matrix(pred)), unname(ref_predict(fit, test[x_names], 2)), tolerance = 1e-6)

  idx <- c(1, 51, 101)
  fit_cls <- lwpls_fit(iris[-idx, 1:4], iris$Species[-idx], num_comp = 2, robust = TRUE)
  raw <- predict(fit_cls, iris[idx, 1:4], type = "raw")[, , 2]
  expect_equal(unname(raw), unname(ref_predict(fit_cls, iris[idx, 1:4], 2)), tolerance = 1e-6)
  # The predicted class indicators still sum to one
  expect_equal(unname(rowSums(raw)), rep(1, 3))
  prob <- predict(fit_cls, iris[idx, 1:4], type = "prob")
  expect_equal(rowSums(prob), rep(1, 3))
})

test_that("robust models in the global limit are partial robust M-regression", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, localization = 1e12, robust = TRUE, scale = FALSE)
  xq <- as.matrix(test[x_names])
  ref <- unname(apply(sweep(xq, 2, fit$x_center), 1, function(q) {
    local_prm(fit$x, fit$y, q, rep(1, nrow(fit$x)), 2)
  })) + fit$y_center
  expect_equal(predict(fit, test[x_names])$.pred, ref, tolerance = 1e-6)
})

test_that("robust models resist outliers in the reference values", {
  n <- 120
  x <- withr::with_seed(11, matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("v", 1:5))))
  truth <- drop(x %*% c(2, -1, 1, 0, 0.5))
  y <- truth + withr::with_seed(12, rnorm(n, sd = 0.2))
  bad <- withr::with_seed(13, sample(100, 20))
  y[bad] <- y[bad] + 20
  rmse <- function(fit) sqrt(mean((predict(fit, x[101:120, ])$.pred - truth[101:120])^2))

  classic <- lwpls_fit(x[1:100, ], y[1:100], num_comp = 3, localization = 2)
  fair <- lwpls_fit(x[1:100, ], y[1:100], num_comp = 3, localization = 2, robust = TRUE)
  hampel <- lwpls_fit(
    x[1:100, ], y[1:100],
    num_comp = 3, localization = 2, robust = TRUE, weight_function = "hampel"
  )
  expect_lt(rmse(fair), rmse(classic) / 3)
  expect_lt(rmse(hampel), rmse(classic) / 3)
})

test_that("robust models are refitted for each number of components", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, robust = TRUE)
  raw <- predict(fit, test[x_names], type = "raw")
  expect_false(anyNA(raw))
  for (a in 1:3) {
    expect_equal(unname(raw[, 1, a]), predict(fit, test[x_names], num_comp = a)$.pred)
  }
})

test_that("robust models handle missing values and degenerate data", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, robust = TRUE)
  new <- test[1:3, x_names]
  new$x1[2] <- NA
  pred <- predict(fit, new)$.pred
  expect_true(is.na(pred[2]))
  expect_equal(pred[-2], predict(fit, new[-2, ])$.pred)

  x <- data.frame(a = rep(1, 10), b = rep(2, 10))
  y <- as.numeric(1:10)
  const <- lwpls_fit(x, y, num_comp = 1, robust = TRUE)
  expect_equal(predict(const, x[1:2, ])$.pred, rep(stats::median(y), 2))

  # One neighbor: no local component can be fitted
  nn <- lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 1, robust = TRUE)
  expect_true(all(is.finite(predict(nn, test[x_names])$.pred)))
})

test_that("edge cases of sparse and robust models", {
  # Orthogonal predictors and y = a: nothing is left to model after the first
  # component, whose weight vector is (1, 0)
  x <- cbind(a = rep(c(-1, 1), 20), b = rep(c(-1, -1, 1, 1), 10))
  lin <- lwpls_fit(x, x[, "a"], num_comp = 2, localization = 1e15, sparsity = 0.1)
  new <- cbind(a = c(0.5, -2), b = c(1, 3))
  raw <- predict(lin, new, type = "raw")
  expect_equal(unname(raw[, 1, 1]), c(0.5, -2))
  expect_equal(unname(raw[, 1, 2]), unname(raw[, 1, 1]))

  # One neighbor: no component, the prediction is that neighbor's outcome
  one <- lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 1, sparsity = 0.2)
  nearest <- lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 1)
  expect_equal(predict(one, test[x_names]), predict(nearest, test[x_names]))

  # Many tied outcomes: more than half of the residuals are zero at first
  tied <- train
  tied$y[1:40] <- stats::median(tied$y)
  fit_tied <- lwpls_fit(tied[x_names], tied$y, num_comp = 2, robust = TRUE)
  expect_true(all(is.finite(predict(fit_tied, test[x_names])$.pred)))
})

test_that("robust covariance-based similarity with duplicated samples", {
  # More than half of the samples at the center of the predictors and outcome
  dup <- train[c(x_names, "y")]
  dup[1:40, ] <- matrix(
    vapply(dup, stats::median, numeric(1)),
    nrow = 40, ncol = ncol(dup), byrow = TRUE
  )
  fit <- lwpls_fit(dup[x_names], dup$y, num_comp = 2, robust = TRUE, similarity = "covariance")
  expect_true(all(is.finite(fit$projection)))
  expect_true(all(is.finite(predict(fit, test[x_names])$.pred)))
  # All samples identical in the predictors: no covariance direction
  const <- data.frame(a = rep(1, 10), b = rep(2, 10))
  y <- as.numeric(1:10)
  fit_const <- lwpls_fit(const, y, num_comp = 1, robust = TRUE, similarity = "covariance")
  expect_equal(unname(fit_const$projection), matrix(0, 2, 1))
  expect_equal(predict(fit_const, const[1:2, ])$.pred, rep(stats::median(y), 2))
})

test_that("robust covariance-based similarity", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, robust = TRUE, similarity = "covariance")
  expect_equal(dim(fit$projection), c(6, 1))
  expect_true(all(is.finite(predict(fit, test[x_names])$.pred)))

  fit_cls <- lwpls_fit(Species ~ ., data = iris, robust = TRUE, similarity = "covariance")
  expect_gt(mean(predict(fit_cls, iris)$.pred_class == iris$Species), 0.9)
})

# Interface -------------------------------------------------------------------

test_that("print method for sparse and robust models", {
  expect_snapshot(lwpls_fit(train[x_names], train$y, sparsity = 0.3, robust = TRUE))
  expect_snapshot(
    lwpls_fit(train[x_names], train$y, robust = TRUE, weight_function = "hampel")
  )
})

test_that("bad sparse and robust settings are rejected", {
  x <- train[x_names]
  y <- train$y
  expect_snapshot(error = TRUE, lwpls_fit(x, y, sparsity = 1))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, sparsity = -0.1))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, robust = "yes"))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, robust = TRUE, weight_function = "huber"))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, robust = TRUE, robust_constant = -1))
  expect_snapshot(
    error = TRUE,
    lwpls_fit(x, y, robust = TRUE, weight_function = "hampel", robust_constant = c(0.9, 0.8, 0.99))
  )
  expect_snapshot(error = TRUE, lwpls_fit(x, y, robust = TRUE, max_iter = 0))

  fit <- lwpls_fit(x, y, robust = TRUE, weight_function = "hampel", robust_constant = c(0.9, 0.95, 0.99))
  expect_equal(fit$robust_constant, c(0.9, 0.95, 0.99))
  expect_equal(lwpls_fit(x, y, robust = TRUE, robust_constant = 2)$robust_constant, 2)
})

test_that("parsnip passes the sparse and robust settings", {
  spec <- lwpls(num_comp = 3, sparsity = 0.3) |>
    parsnip::set_mode("regression") |>
    parsnip::set_engine("lwpls", robust = TRUE, weight_function = "hampel")
  expect_snapshot(parsnip::translate(spec))

  fit <- parsnip::fit_xy(spec, x = train[x_names], y = train$y)
  expect_true(fit$fit$robust)
  expect_equal(fit$fit$sparsity, 0.3)
  direct <- lwpls_fit(
    train[x_names], train$y,
    num_comp = 3, sparsity = 0.3, robust = TRUE, weight_function = "hampel"
  )
  expect_equal(predict(fit, test[x_names]), predict(direct, test[x_names]))

  mp <- parsnip::multi_predict(fit, test[1:4, x_names], num_comp = 1:3)
  for (a in 1:3) {
    by_comp <- vapply(mp$.pred, function(x) x$.pred[x$num_comp == a], numeric(1))
    expect_equal(by_comp, predict(direct, test[1:4, x_names], num_comp = a)$.pred)
  }

  expect_snapshot(
    error = TRUE,
    lwpls(sparsity = 2) |>
      parsnip::set_mode("regression") |>
      parsnip::fit(y ~ ., data = train[c(x_names, "y")])
  )
})

test_that("sparsity parameter", {
  param <- sparsity()
  expect_s3_class(param, "quant_param")
  expect_equal(param$range, list(lower = 0, upper = 0.9))
  expect_snapshot(sparsity())
})

test_that("tune_grid() with sparse robust models", {
  skip_on_cran()
  skip_if_not_installed("tune")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")

  folds <- withr::with_seed(1, rsample::vfold_cv(train[c(x_names, "y")], v = 3))
  spec <- lwpls(num_comp = tune::tune(), sparsity = tune::tune()) |>
    parsnip::set_mode("regression") |>
    parsnip::set_engine("lwpls", robust = TRUE)
  res <- tune::tune_grid(
    workflows::workflow(y ~ ., spec),
    folds,
    grid = expand.grid(num_comp = 1:2, sparsity = c(0, 0.4)),
    metrics = yardstick::metric_set(trimmed_rmse, yardstick::rmse)
  )
  perf <- tune::collect_metrics(res)
  expect_equal(nrow(perf), 8)
  expect_true(all(is.finite(perf$mean)))
  expect_setequal(perf$.metric, c("trimmed_rmse", "rmse"))
})

# Trimmed RMSE ----------------------------------------------------------------

test_that("trimmed RMSE", {
  df <- data.frame(truth = c(1, 2, 3, 4, 5), estimate = c(1.1, 2.1, 2.9, 4.2, 15))
  res <- trimmed_rmse(df, truth, estimate)
  expect_named(res, c(".metric", ".estimator", ".estimate"))
  expect_equal(res$.metric, "trimmed_rmse")
  expect_equal(res$.estimate, sqrt(mean(c(0.1, 0.1, 0.1, 0.2)^2)))
  expect_equal(
    trimmed_rmse_vec(df$truth, df$estimate, trim = 0),
    yardstick::rmse_vec(df$truth, df$estimate)
  )
  # Case weights: a weight of two counts an observation twice
  expect_equal(
    trimmed_rmse_vec(c(1, 2, 3), c(1.5, 2, 3), trim = 0, case_weights = c(2, 1, 1)),
    sqrt(2 * 0.25 / 4)
  )
  # At least one observation is kept
  expect_equal(trimmed_rmse_vec(c(1, 2), c(1.5, 3), trim = 0.9), 0.5)

  na <- c(1, NA, 3)
  expect_equal(trimmed_rmse_vec(na, c(1.5, 2, 3), trim = 0), sqrt(0.25 / 2))
  expect_true(is.na(trimmed_rmse_vec(na, c(1.5, 2, 3), na_rm = FALSE)))
  expect_true(is.na(trimmed_rmse_vec(c(NA_real_, NA_real_), c(1, 2))))
  expect_snapshot(error = TRUE, trimmed_rmse_vec(1:3, 1:3, trim = 1))

  grouped <- data.frame(g = c("a", "a", "b", "b"), truth = 1:4, estimate = c(1, 3, 3, 4))
  by_group <- trimmed_rmse(dplyr::group_by(grouped, g), truth, estimate, trim = 0)
  expect_equal(by_group$.estimate, c(sqrt(0.5), 0))
})
