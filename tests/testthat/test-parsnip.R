dat <- sim_reg()
x_names <- paste0("x", 1:6)
train <- dat[1:90, c(x_names, "y", "y2")]
test <- dat[91:120, c(x_names, "y", "y2")]

# Model specification ---------------------------------------------------------

test_that("model specification", {
  spec <- lwpls(num_comp = 3, localization = 0.5)
  expect_s3_class(spec, c("lwpls", "model_spec"))
  expect_equal(spec$engine, "lwpls")
  expect_snapshot(spec)
  expect_snapshot(
    lwpls(num_comp = 3, neighbors = 50) |>
      parsnip::set_mode("classification") |>
      parsnip::set_engine("lwpls", scale = FALSE) |>
      parsnip::translate()
  )
})

test_that("update method", {
  spec <- lwpls(num_comp = 3)
  expect_snapshot(update(spec, localization = 0.25))
  expect_snapshot(update(spec, localization = 0.25, fresh = TRUE))

  param <- tibble::tibble(num_comp = 5, localization = 2)
  upd <- update(spec, param)
  expect_equal(rlang::eval_tidy(upd$args$num_comp), 5)
  expect_equal(rlang::eval_tidy(upd$args$localization), 2)

  expect_snapshot(error = TRUE, update(spec, param, scale = FALSE))
})

test_that("argument checks happen at fit time", {
  expect_snapshot(
    error = TRUE,
    lwpls(localization = -1) |>
      parsnip::set_mode("regression") |>
      parsnip::fit(y ~ ., data = train)
  )
  expect_snapshot(
    error = TRUE,
    lwpls(num_comp = 0) |>
      parsnip::set_mode("regression") |>
      parsnip::fit(y ~ ., data = train)
  )
})

test_that("engine dependencies", {
  spec <- lwpls() |> parsnip::set_mode("regression")
  expect_equal(parsnip::required_pkgs(spec), c("parsnip", "lwpls"))
})

# Regression ------------------------------------------------------------------

test_that("regression fit and predictions", {
  spec <- lwpls(num_comp = 3, localization = 0.5) |> parsnip::set_mode("regression")

  fit_form <- parsnip::fit(spec, y ~ . - y2, data = train)
  fit_xy <- parsnip::fit_xy(spec, x = train[x_names], y = train$y)
  expect_s3_class(fit_form, "_lwpls_fit")
  expect_s3_class(fit_form$fit, "lwpls_fit")

  direct <- lwpls_fit(train[x_names], train$y, num_comp = 3, localization = 0.5)
  expected <- predict(direct, test[x_names])

  pred_form <- predict(fit_form, test)
  expect_named(pred_form, ".pred")
  expect_equal(pred_form, expected)
  expect_equal(predict(fit_xy, test[x_names]), expected)
  expect_equal(nrow(parsnip::augment(fit_form, test)), nrow(test))

  raw <- predict(fit_xy, test[x_names], type = "raw")
  expect_equal(dim(raw), c(nrow(test), 1, 3))
})

test_that("the similarity index is passed on", {
  spec <- lwpls(num_comp = 2, similarity = "covariance") |>
    parsnip::set_mode("regression")
  expect_snapshot(parsnip::translate(spec))
  fit <- parsnip::fit_xy(spec, x = train[x_names], y = train$y)
  expect_equal(fit$fit$similarity, "covariance")
  direct <- lwpls_fit(train[x_names], train$y, num_comp = 2, similarity = "covariance")
  expect_equal(predict(fit, test[x_names]), predict(direct, test[x_names]))

  expect_equal(
    rlang::eval_tidy(update(spec, similarity = "euclidean")$args$similarity),
    "euclidean"
  )
  expect_snapshot(
    error = TRUE,
    lwpls(similarity = "cosine") |>
      parsnip::set_mode("regression") |>
      parsnip::fit(y ~ ., data = train)
  )
})

test_that("engine arguments are passed on", {
  spec <- lwpls(num_comp = 2) |>
    parsnip::set_mode("regression") |>
    parsnip::set_engine("lwpls", scale = FALSE)
  fit <- parsnip::fit_xy(spec, x = train[x_names], y = train$y)
  expect_false(fit$fit$scale)
})

test_that("multivariate regression", {
  spec <- lwpls(num_comp = 3) |> parsnip::set_mode("regression")
  fit <- parsnip::fit(spec, cbind(y, y2) ~ ., data = train)
  pred <- predict(fit, test)
  expect_named(pred, c(".pred_y", ".pred_y2"))

  direct <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 3)
  expect_equal(pred, predict(direct, test[x_names]))

  mp <- parsnip::multi_predict(fit, test[1:3, ], num_comp = 1:3)
  expect_named(mp$.pred[[1]], c("num_comp", ".pred_y", ".pred_y2"))
})

test_that("multi_predict() for regression", {
  spec <- lwpls(num_comp = 4) |> parsnip::set_mode("regression")
  fit <- parsnip::fit(spec, y ~ . - y2, data = train)
  expect_true(parsnip::has_multi_predict(fit))
  expect_equal(parsnip::multi_predict_args(fit), "num_comp")

  mp <- parsnip::multi_predict(fit, test, num_comp = c(3, 1, 2))
  expect_s3_class(mp, "tbl_df")
  expect_named(mp, ".pred")
  expect_equal(nrow(mp), nrow(test))
  expect_true(all(vapply(mp$.pred, nrow, integer(1)) == 3))
  expect_named(mp$.pred[[1]], c("num_comp", ".pred"))
  expect_equal(mp$.pred[[1]]$num_comp, 1:3)

  for (a in 1:3) {
    by_comp <- vapply(mp$.pred, function(x) x$.pred[x$num_comp == a], numeric(1))
    expect_equal(by_comp, predict(fit$fit, test, num_comp = a)$.pred)
  }

  # Default: the fitted number of components
  mp_default <- parsnip::multi_predict(fit, test[1:2, ])
  expect_equal(mp_default$.pred[[1]]$num_comp, 4L)

  expect_snapshot(error = TRUE, parsnip::multi_predict(fit, newdata = test))
  expect_snapshot(error = TRUE, parsnip::multi_predict(fit, test, type = "prob"))
  expect_snapshot(error = TRUE, parsnip::multi_predict(fit, test, num_comp = -1))
})

# Classification --------------------------------------------------------------

test_that("classification fit and predictions", {
  spec <- lwpls(num_comp = 2, localization = 0.5) |> parsnip::set_mode("classification")
  fit <- parsnip::fit(spec, Species ~ ., data = iris)
  new <- iris[c(1, 51, 101), ]

  cls <- predict(fit, new)
  prob <- predict(fit, new, type = "prob")
  expect_named(cls, ".pred_class")
  expect_s3_class(cls$.pred_class, "factor")
  expect_equal(levels(cls$.pred_class), levels(iris$Species))
  expect_named(prob, paste0(".pred_", levels(iris$Species)))

  direct <- lwpls_fit(iris[1:4], iris$Species, num_comp = 2, localization = 0.5)
  expect_equal(prob, predict(direct, new[1:4], type = "prob"))
  expect_equal(cls, predict(direct, new[1:4]))
  expect_snapshot(error = TRUE, predict(fit, new, type = "numeric"))
})

test_that("multi_predict() for classification", {
  spec <- lwpls(num_comp = 3) |> parsnip::set_mode("classification")
  fit <- parsnip::fit(spec, Species ~ ., data = iris)
  new <- iris[c(1, 51, 101), ]

  mp_cls <- parsnip::multi_predict(fit, new, num_comp = 1:3)
  expect_named(mp_cls$.pred[[1]], c("num_comp", ".pred_class"))

  mp_prob <- parsnip::multi_predict(fit, new, type = "prob", num_comp = 1:3)
  expect_named(
    mp_prob$.pred[[1]],
    c("num_comp", paste0(".pred_", levels(iris$Species)))
  )
  expect_equal(
    mp_prob$.pred[[3]][2, -1],
    predict(fit$fit, new[3, ], type = "prob", num_comp = 2)
  )
})
