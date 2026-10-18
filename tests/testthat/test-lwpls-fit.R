dat <- sim_reg()
x_names <- paste0("x", 1:6)
train <- dat[1:90, ]
test <- dat[91:120, ]

# Algorithm -------------------------------------------------------------------

test_that("predictions match Kaneko's reference implementation", {
  for (lambda in c(0.25, 1, 4)) {
    fit <- lwpls_fit(train[x_names], train$y, num_comp = 4, localization = lambda)
    raw <- predict(fit, test[x_names], type = "raw")
    ref <- kaneko_predict(as.matrix(train[x_names]), train$y, as.matrix(test[x_names]), 4, lambda)
    expect_equal(unname(raw[, 1, ]), ref, tolerance = 1e-10)
  }
})

test_that("predictions match the reference without scaling", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, localization = 0.5, scale = FALSE)
  raw <- predict(fit, test[x_names], type = "raw")
  ref <- kaneko_predict(
    as.matrix(train[x_names]), train$y, as.matrix(test[x_names]), 3, 0.5,
    scale = FALSE
  )
  expect_equal(unname(raw[, 1, ]), ref, tolerance = 1e-10)
})

test_that("covariance-based similarity matches the CbLW-PLS reference", {
  for (lambda in c(0.25, 1, 4)) {
    fit <- lwpls_fit(
      train[x_names], train$y,
      num_comp = 4, localization = lambda, similarity = "covariance"
    )
    raw <- predict(fit, test[x_names], type = "raw")
    ref <- kaneko_predict(
      as.matrix(train[x_names]), train$y, as.matrix(test[x_names]), 4, lambda,
      covariance = TRUE
    )
    expect_equal(unname(raw[, 1, ]), ref, tolerance = 1e-10)
  }

  euclidean <- lwpls_fit(train[x_names], train$y, num_comp = 4)
  expect_null(euclidean$projection)
  expect_false(isTRUE(all.equal(
    predict(fit, test[x_names]),
    predict(euclidean, test[x_names])
  )))
})

test_that("covariance-based similarity with several outcomes and classes", {
  y <- as.matrix(train[c("y", "y2")])
  fit <- lwpls_fit(train[x_names], y, num_comp = 3, similarity = "covariance")
  xy <- crossprod(fit$x, fit$y)
  expect_equal(fit$projection, xy / sqrt(sum(xy^2)))
  expect_equal(dim(fit$projection), c(6, 2))
  pred <- predict(fit, test[x_names])
  expect_named(pred, c(".pred_y", ".pred_y2"))
  expect_true(all(is.finite(as.matrix(pred))))

  fit_cls <- lwpls_fit(Species ~ ., data = iris, num_comp = 2, similarity = "covariance")
  prob <- predict(fit_cls, iris[seq(1, 150, by = 10), ], type = "prob")
  expect_equal(rowSums(prob), rep(1, nrow(prob)))
  expect_gt(mean(predict(fit_cls, iris)$.pred_class == iris$Species), 0.9)
})

test_that("covariance-based similarity without covariance gives global PLS weights", {
  x <- data.frame(a = c(1, -1, 1, -1), b = c(1, 1, -1, -1))
  y <- c(1, 2, 2, 1)
  fit <- lwpls_fit(x, y, num_comp = 1, similarity = "covariance")
  expect_equal(unname(fit$projection), matrix(0, 2, 1))
  expect_equal(predict(fit, x)$.pred, rep(mean(y), 4))
})

test_that("`neighbors` restricts each local model to the nearest samples", {
  k <- 30
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, localization = 1, neighbors = k)
  pred <- predict(fit, test[x_names])$.pred

  mu <- colMeans(train[x_names])
  s <- apply(train[x_names], 2, sd)
  xs <- scale(as.matrix(train[x_names]), mu, s)
  ys <- (train$y - mean(train$y)) / sd(train$y)
  qs <- scale(as.matrix(test[x_names]), mu, s)
  ref <- vapply(seq_len(nrow(qs)), function(i) {
    d <- sqrt(colSums((t(xs) - qs[i, ])^2))
    nn <- order(d)[1:k]
    kaneko_lwpls(xs[nn, ], ys[nn], qs[i, , drop = FALSE], 3, 1)[1, 3]
  }, numeric(1))
  expect_equal(pred, ref * sd(train$y) + mean(train$y), tolerance = 1e-10)
})

test_that("one neighbor gives nearest-neighbor predictions", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 1)
  pred <- predict(fit, test[x_names])$.pred
  xs <- scale(as.matrix(train[x_names]))
  qs <- scale(
    as.matrix(test[x_names]),
    attr(xs, "scaled:center"),
    attr(xs, "scaled:scale")
  )
  nn <- apply(qs, 1, function(q) which.min(colSums((t(xs) - q)^2)))
  expect_equal(pred, train$y[nn])
})

test_that("a very large localization gives global PLS", {
  skip_if_not_installed("pls")
  y <- as.matrix(train[c("y", "y2")])
  fit <- lwpls_fit(train[x_names], y, num_comp = 4, localization = 1e12)
  raw <- predict(fit, test[x_names], type = "raw")

  xs <- scale(as.matrix(train[x_names]))
  ys <- scale(y)
  qs <- scale(
    as.matrix(test[x_names]),
    attr(xs, "scaled:center"),
    attr(xs, "scaled:scale")
  )
  pls_fit <- pls::plsr(ys ~ xs, ncomp = 4, method = "kernelpls")
  ref <- predict(pls_fit, newdata = list(xs = qs), ncomp = 1:4)
  ref <- sweep(ref, 2, attr(ys, "scaled:scale"), "*")
  ref <- sweep(ref, 2, attr(ys, "scaled:center"), "+")
  expect_equal(unname(raw), unname(ref), tolerance = 1e-7)
})

test_that("predictions are stable for very local models", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 5, localization = 2^-9)
  raw <- predict(fit, test[x_names], type = "raw")
  expect_true(all(is.finite(raw)))
  expect_true(all(raw >= min(train$y) - 1 & raw <= max(train$y) + 1))
})

test_that("degenerate local models fall back to the weighted mean", {
  x <- data.frame(a = rep(1, 10), b = rep(2, 10))
  y <- as.numeric(1:10)
  expect_no_warning(fit <- lwpls_fit(x, y, num_comp = 1))
  expect_equal(predict(fit, x[1:2, ])$.pred, rep(mean(y), 2))
})

test_that("wide data (p > n) is supported", {
  wide <- withr::with_seed(2, matrix(rnorm(20 * 50), 20, 50))
  colnames(wide) <- paste0("v", 1:50)
  y <- wide[, 1] + rnorm(20, sd = 0.1)
  expect_snapshot(fit <- lwpls_fit(wide[1:15, ], y[1:15], num_comp = 20))
  expect_equal(fit$num_comp, 14L)
  raw <- predict(fit, wide[16:20, ], type = "raw")
  expect_equal(dim(raw), c(5, 1, 14))
  expect_true(all(is.finite(raw)))
})

test_that("constant predictors are handled", {
  x <- train[x_names]
  x$const <- 3
  fit <- lwpls_fit(x, train$y, num_comp = 3)
  new <- test[x_names]
  new$const <- 3
  ref <- lwpls_fit(train[x_names], train$y, num_comp = 3)
  expect_equal(predict(fit, new), predict(ref, test[x_names]))
})

# Interfaces ------------------------------------------------------------------

test_that("all interfaces give the same model", {
  fit_df <- lwpls_fit(train[x_names], train$y, num_comp = 3)
  fit_mat <- lwpls_fit(as.matrix(train[x_names]), train$y, num_comp = 3)
  fit_form <- lwpls_fit(y ~ ., data = train[c(x_names, "y")], num_comp = 3)

  expected <- predict(fit_df, test[x_names])
  expect_equal(predict(fit_mat, as.matrix(test[x_names])), expected)
  expect_equal(predict(fit_form, test), expected)

  skip_if_not_installed("recipes")
  rec <- recipes::recipe(y ~ ., data = train[c(x_names, "y")])
  fit_rec <- lwpls_fit(rec, data = train[c(x_names, "y")], num_comp = 3)
  expect_equal(predict(fit_rec, test), expected)
})

test_that("the formula interface creates dummy variables", {
  df <- train[c(x_names, "y")]
  df$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(df)))
  fit <- lwpls_fit(y ~ ., data = df, num_comp = 2)
  expect_equal(ncol(fit$x), 8)
  new <- test
  new$grp <- factor(rep("b", nrow(test)), levels = c("a", "b", "c"))
  expect_equal(nrow(predict(fit, new)), nrow(test))
})

test_that("multiple numeric outcomes are supported", {
  fit <- lwpls_fit(y + y2 ~ ., data = train, num_comp = 3)
  pred <- predict(fit, test)
  expect_named(pred, c(".pred_y", ".pred_y2"))
  expect_equal(nrow(pred), nrow(test))

  fit_xy <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 3)
  expect_equal(predict(fit_xy, test[x_names]), pred)
})

# Classification --------------------------------------------------------------

test_that("classification predictions are consistent", {
  fit <- lwpls_fit(Species ~ ., data = iris, num_comp = 2, localization = 0.5)
  expect_equal(fit$mode, "classification")

  new <- iris[seq(1, 150, by = 10), ]
  cls <- predict(fit, new)
  prob <- predict(fit, new, type = "prob")
  raw <- predict(fit, new, type = "raw")

  expect_named(cls, ".pred_class")
  expect_equal(levels(cls$.pred_class), levels(iris$Species))
  expect_named(prob, paste0(".pred_", levels(iris$Species)))
  expect_equal(rowSums(prob), rep(1, nrow(new)))
  expect_true(all(prob >= 0 & prob <= 1))
  # Predicted class indicators sum to one
  expect_equal(apply(raw, c(1, 3), sum), matrix(1, nrow(new), 2, dimnames = list(NULL, 1:2)))
  expect_equal(
    as.character(cls$.pred_class),
    levels(iris$Species)[max.col(as.matrix(prob), ties.method = "first")]
  )
  expect_gt(mean(cls$.pred_class == new$Species), 0.9)
})

test_that("classification with two classes works", {
  df <- iris[iris$Species != "setosa", ]
  df$Species <- droplevels(df$Species)
  fit <- lwpls_fit(df[1:4], df$Species, num_comp = 2)
  prob <- predict(fit, df[1:5, 1:4], type = "prob")
  expect_named(prob, c(".pred_versicolor", ".pred_virginica"))
  expect_equal(rowSums(prob), rep(1, 5))
})

# Prediction options ----------------------------------------------------------

test_that("`num_comp` can be changed at prediction time", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2)
  raw <- predict(fit, test[x_names], type = "raw", num_comp = 5)
  expect_equal(dim(raw), c(nrow(test), 1, 5))
  expect_equal(predict(fit, test[x_names])$.pred, unname(raw[, 1, 2]))
  expect_equal(predict(fit, test[x_names], num_comp = 4)$.pred, unname(raw[, 1, 4]))

  # Beyond the maximum (p = 6), predictions stay at the maximum
  expect_equal(
    predict(fit, test[x_names], num_comp = 9),
    predict(fit, test[x_names], num_comp = 6)
  )
})

test_that("missing values in new data give missing predictions", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2)
  new <- test[1:4, x_names]
  new$x2[2] <- NA
  pred <- predict(fit, new)
  expect_equal(nrow(pred), 4)
  expect_true(is.na(pred$.pred[2]))
  expect_equal(pred$.pred[-2], predict(fit, new[-2, ])$.pred)

  all_missing <- new
  all_missing$x1 <- NA_real_
  expect_equal(predict(fit, all_missing)$.pred, rep(NA_real_, 4))

  fit_cls <- lwpls_fit(iris[1:4], iris$Species)
  new_cls <- iris[c(1, 51), 1:4]
  new_cls$Sepal.Width[1] <- NA
  expect_true(is.na(predict(fit_cls, new_cls)$.pred_class[1]))
  expect_true(all(is.na(predict(fit_cls, new_cls, type = "prob")[1, ])))
})

test_that("print method", {
  expect_snapshot(lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 50))
  expect_snapshot(lwpls_fit(Species ~ ., data = iris, scale = FALSE))
  expect_snapshot(lwpls_fit(y + y2 ~ ., data = train, num_comp = 2))
  expect_snapshot(lwpls_fit(y ~ ., data = train, similarity = "covariance"))
})

# Input validation ------------------------------------------------------------

test_that("bad hyperparameters are rejected", {
  x <- train[x_names]
  y <- train$y
  expect_snapshot(error = TRUE, lwpls_fit(x, y, num_comp = 0))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, num_comp = 1.5))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, num_comp = 1:2))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, localization = -1))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, localization = Inf))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, localization = "small"))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, neighbors = 0))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, similarity = "mahalanobis"))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, scale = "yes"))
  expect_snapshot(error = TRUE, lwpls_fit(x, y, num_comps = 2))
  expect_snapshot(fit <- lwpls_fit(x, y, neighbors = 500))
  expect_equal(fit$neighbors, nrow(x))
})

test_that("bad data are rejected", {
  x <- train[x_names]
  y <- train$y
  expect_snapshot(error = TRUE, lwpls_fit(list(1)))
  expect_snapshot(error = TRUE, lwpls_fit(matrix(rnorm(20), 10), y[1:10]))
  expect_snapshot(error = TRUE, lwpls_fit(iris[5:4], iris$Sepal.Length))
  expect_snapshot(error = TRUE, lwpls_fit(x[1, ], y[1]))
  expect_snapshot(error = TRUE, lwpls_fit(x[0], y))

  x_na <- x
  x_na$x1[3] <- NA
  expect_snapshot(error = TRUE, lwpls_fit(x_na, y))
  y_na <- y
  y_na[3] <- NA
  expect_snapshot(error = TRUE, lwpls_fit(x, y_na))
  expect_snapshot(error = TRUE, lwpls_fit(x, factor(ifelse(y_na > 1, "a", "b"))))
  expect_snapshot(error = TRUE, lwpls_fit(x, as.character(y > 1)))
  expect_snapshot(error = TRUE, lwpls_fit(x, factor(rep("a", nrow(x)))))
  expect_snapshot(
    error = TRUE,
    lwpls_fit(x, data.frame(a = y, b = factor(y > 1)))
  )
})

test_that("bad prediction arguments are rejected", {
  fit <- lwpls_fit(train[x_names], train$y)
  expect_snapshot(error = TRUE, predict(fit, test[x_names], type = "class"))
  expect_snapshot(error = TRUE, predict(fit, test[x_names], num_comp = 0))
  expect_snapshot(error = TRUE, predict(fit, test[x_names], num_comp = 1:2))
  expect_snapshot(error = TRUE, predict(fit, test[x_names], foo = 1))
  expect_error(predict(fit, test[x_names[-1]]))
})
