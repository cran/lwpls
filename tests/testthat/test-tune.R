test_that("localization parameter", {
  param <- localization()
  expect_s3_class(param, "quant_param")
  expect_equal(param$range, list(lower = -9, upper = 5))
  expect_equal(dials::value_seq(param, 3, original = TRUE), 2^c(-9, -2, 5))
  expect_snapshot(localization())
})

test_that("similarity parameter", {
  param <- similarity()
  expect_s3_class(param, "qual_param")
  expect_equal(param$values, c("euclidean", "covariance"))
  expect_equal(similarity("covariance")$values, "covariance")
  expect_snapshot(similarity())
})

test_that("tunable parameters", {
  skip_if_not_installed("tune")
  spec <- lwpls(
    num_comp = tune::tune(),
    localization = tune::tune(),
    neighbors = tune::tune(),
    similarity = tune::tune(),
    sparsity = tune::tune()
  ) |>
    parsnip::set_mode("regression")

  tunable <- generics::tunable(spec)
  expect_equal(tunable$name, c("num_comp", "localization", "neighbors", "similarity", "sparsity"))

  params <- hardhat::extract_parameter_set_dials(spec)
  expect_equal(params$id, c("num_comp", "localization", "neighbors", "similarity", "sparsity"))
  expect_equal(params$object[[2]]$range, localization()$range)
  expect_equal(params$object[[3]]$range, list(lower = 10L, upper = 200L))
  expect_equal(params$object[[4]]$values, values_similarity)
  expect_equal(params$object[[5]]$range, sparsity()$range)
})

test_that("min_grid() uses the submodel trick for num_comp", {
  skip_if_not_installed("tune")
  spec <- lwpls(num_comp = tune::tune(), localization = tune::tune()) |>
    parsnip::set_mode("regression")
  grid <- expand.grid(num_comp = 1:4, localization = c(0.5, 1))
  res <- min_grid(spec, grid)
  expect_equal(nrow(res), 2)
  expect_equal(res$num_comp, c(4L, 4L))
  expect_equal(res$.submodels[[1]], list(num_comp = 1:3))
})

test_that("tune_grid() works and matches direct resampling", {
  skip_on_cran()
  skip_if_not_installed("tune")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")

  dat <- sim_reg(n = 90)[c(paste0("x", 1:6), "y")]
  folds <- withr::with_seed(1, rsample::vfold_cv(dat, v = 3))
  metrics <- yardstick::metric_set(yardstick::rmse)

  spec <- lwpls(num_comp = tune::tune(), localization = tune::tune()) |>
    parsnip::set_mode("regression")
  wflow <- workflows::workflow(y ~ ., spec)
  grid <- expand.grid(num_comp = 1:4, localization = c(0.5, 2))

  res <- tune::tune_grid(wflow, folds, grid = grid, metrics = metrics)
  perf <- tune::collect_metrics(res)
  expect_equal(nrow(perf), nrow(grid))
  expect_true(all(is.finite(perf$mean)))

  direct <- tune::fit_resamples(
    workflows::workflow(y ~ ., lwpls(num_comp = 2, localization = 0.5) |> parsnip::set_mode("regression")),
    folds,
    metrics = metrics
  )
  expect_equal(
    perf$mean[perf$num_comp == 2 & perf$localization == 0.5],
    tune::collect_metrics(direct)$mean
  )
})

test_that("tune_grid() can compare similarity indexes", {
  skip_on_cran()
  skip_if_not_installed("tune")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")

  dat <- sim_reg(n = 90)[c(paste0("x", 1:6), "y")]
  folds <- withr::with_seed(1, rsample::vfold_cv(dat, v = 3))
  spec <- lwpls(num_comp = tune::tune(), similarity = tune::tune()) |>
    parsnip::set_mode("regression")
  grid <- expand.grid(num_comp = 1:3, similarity = values_similarity, stringsAsFactors = FALSE)

  res <- tune::tune_grid(
    workflows::workflow(y ~ ., spec),
    folds,
    grid = grid,
    metrics = yardstick::metric_set(yardstick::rmse)
  )
  perf <- tune::collect_metrics(res)
  expect_equal(nrow(perf), 6)
  expect_setequal(perf$similarity, values_similarity)
  expect_true(all(is.finite(perf$mean)))
})

test_that("tune_grid() works for classification", {
  skip_on_cran()
  skip_if_not_installed("tune")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")

  folds <- withr::with_seed(1, rsample::vfold_cv(iris, v = 3, strata = Species))
  spec <- lwpls(num_comp = tune::tune()) |> parsnip::set_mode("classification")
  res <- tune::tune_grid(
    workflows::workflow(Species ~ ., spec),
    folds,
    grid = data.frame(num_comp = 1:3),
    metrics = yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)
  )
  perf <- tune::collect_metrics(res)
  expect_equal(nrow(perf), 6)
  expect_true(all(perf$mean > 0.8))
})
