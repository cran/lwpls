dat <- sim_reg(n = 80, seed = 3)
x_names <- paste0("x", 1:6)
train <- dat[1:65, ]
test <- dat[66:80, ]

# Diagnostics -----------------------------------------------------------------

test_that("local models in the global limit are those of PLS", {
  skip_if_not_installed("pls")
  a <- 3
  fit <- lwpls_fit(train[x_names], train$y, num_comp = a, localization = 1e12)
  loc <- lwpls_local(fit, test[x_names])

  x <- as.matrix(train[x_names])
  pls_fit <- pls::plsr(train$y ~ x, ncomp = a, scale = TRUE, method = "oscorespls")
  new_x <- as.matrix(test[x_names])

  # Regression coefficients on the original scale (pls gives them for the
  # scaled predictors)
  coefs <- loc$coefficients[loc$coefficients$.row == 1, ]
  expect_equal(
    coefs$coefficient,
    unname(drop(coef(pls_fit, ncomp = a)) / pls_fit$scale),
    tolerance = 1e-8
  )
  expect_true(all(coefs$selected))

  # VIP
  W <- pls_fit$loading.weights
  ssy <- drop(pls_fit$Yloadings)^2 * colSums(pls_fit$scores^2)
  vip <- sqrt(ncol(x) * drop(W^2 %*% ssy) / sum(ssy))
  expect_equal(coefs$vip, unname(vip), tolerance = 1e-8)

  # T2 and Q of the queries, and their limits
  n <- nrow(x)
  # pls scales the predictors first; Xmeans are the means of the scaled data
  scaled <- function(z) sweep(sweep(z, 2, pls_fit$scale, "/"), 2, pls_fit$Xmeans)
  scores <- unclass(pls_fit$scores)
  new_scores <- scaled(new_x) %*% pls_fit$projection
  score_var <- apply(scores, 2, stats::var)
  t2 <- rowSums(sweep(new_scores^2, 2, score_var, "/"))
  expect_equal(loc$reliability$t2, unname(t2), tolerance = 1e-6)
  expect_equal(
    loc$reliability$t2_limit_theoretical,
    rep(a * (n^2 - 1) / (n * (n - a)) * stats::qf(0.95, a, n - a), nrow(test)),
    tolerance = 1e-6
  )

  P <- unclass(pls_fit$loadings)
  q_new <- rowSums((scaled(new_x) - new_scores %*% t(P))^2)
  expect_equal(loc$reliability$q, unname(q_new), tolerance = 1e-6)
  q_train <- rowSums((scaled(x) - scores %*% t(P))^2)
  q_mean <- mean(q_train)
  q_var <- mean((q_train - q_mean)^2)
  q_limit <- q_var / (2 * q_mean) * stats::qchisq(0.95, 2 * q_mean^2 / q_var)
  expect_equal(loc$reliability$q_limit_theoretical, rep(q_limit, nrow(test)), tolerance = 1e-6)

  # Empirical limits and percentile ranks among the training samples
  t2_train <- rowSums(sweep(scores^2, 2, score_var, "/"))
  expect_equal(
    loc$reliability$t2_limit_empirical,
    rep(unname(stats::quantile(t2_train, 0.95, type = 1)), nrow(test)),
    tolerance = 1e-6
  )
  expect_equal(
    loc$reliability$q_limit_empirical,
    rep(unname(stats::quantile(q_train, 0.95, type = 1)), nrow(test)),
    tolerance = 1e-6
  )
  expect_equal(loc$reliability$t2_percentile, vapply(t2, function(v) mean(t2_train <= v), numeric(1)), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(loc$reliability$q_percentile, vapply(q_new, function(v) mean(q_train <= v), numeric(1)), tolerance = 1e-6, ignore_attr = TRUE)

  # The default ratios use the empirical limits
  expect_equal(loc$reliability$t2_ratio, loc$reliability$t2 / loc$reliability$t2_limit_empirical)
  theory <- lwpls_local(fit, test[x_names], limits = "theoretical")
  expect_equal(theory$reliability$q_ratio, theory$reliability$q / theory$reliability$q_limit_theoretical)
  strict <- lwpls_local(fit, test[x_names], limits = "theoretical", level = 0.99)
  expect_true(all(strict$reliability$t2_limit_theoretical > theory$reliability$t2_limit_theoretical))
  expect_true(all(strict$reliability$t2_limit_empirical >= loc$reliability$t2_limit_empirical))

  expect_equal(loc$reliability$n_eff, rep(n, nrow(test)), tolerance = 1e-6)
  expect_equal(loc$reliability$.pred, predict(fit, test[x_names])$.pred)
})

test_that("reliability flags extrapolations", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 3, localization = 1)
  far <- test[1:3, x_names] + 6
  loc <- lwpls_local(fit, rbind(test[x_names], far))
  rel <- loc$reliability
  expect_true(all(rel$status[16:18] != "inside"))
  expect_gt(mean(rel$status[1:15] == "inside"), 0.7)
  expect_equal(rel$t2_ratio, rel$t2 / rel$t2_limit_empirical)
  expect_true(all(rel$nearest_distance[16:18] > max(rel$nearest_distance[1:15])))

  # Too few neighbors for the T2 limit
  few <- lwpls_fit(train[x_names], train$y, num_comp = 3, neighbors = 4)
  rel_few <- lwpls_local(few, test[x_names])$reliability
  expect_true(all(rel_few$status == "few neighbors"))
  expect_true(all(is.na(rel_few$t2_limit_theoretical) & !is.nan(rel_few$t2_limit_theoretical)))
  expect_true(all(is.na(rel_few$t2_ratio) & !is.nan(rel_few$t2_ratio)))
  expect_true(all(is.na(rel_few$t2_percentile) & is.na(rel_few$q_limit_empirical)))

  # A single neighbor: no local component, the prediction is that neighbor's
  one <- lwpls_fit(train[x_names], train$y, num_comp = 3, neighbors = 1)
  loc_one <- lwpls_local(one, test[x_names])
  expect_equal(loc_one$reliability$num_comp, rep(0L, nrow(test)))
  expect_equal(loc_one$reliability$.pred, predict(one, test[x_names])$.pred)
  expect_true(all(loc_one$coefficients$coefficient == 0))
  expect_equal(nrow(loc_one$selection), 0)
})

test_that("sparse models: selected predictors and coefficients agree", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, sparsity = 0.6)
  loc <- lwpls_local(fit, test[x_names])
  coefs <- loc$coefficients
  expect_true(all(coefs$coefficient[!coefs$selected] == 0))
  expect_true(any(!coefs$selected))
  expect_true(all(coefs$vip[!coefs$selected] == 0))

  sel <- loc$selection
  expect_named(sel, c("component", "predictor", "position", "frequency", "queries"))
  expect_true(all(sel$frequency >= 0 & sel$frequency <= 1))
  expect_true(all(sel$queries <= nrow(test)))
  # A predictor is selected by a model if any of its components selects it
  first <- tapply(coefs$selected, coefs$predictor, mean)
  expect_true(all(first >= sel$frequency[sel$component == 1] - 1e-12))
})

test_that("robust weights point to the outliers", {
  bad <- train
  bad$y[1:6] <- bad$y[1:6] + 8
  fit <- lwpls_fit(bad[x_names], bad$y, num_comp = 3, localization = 2, robust = TRUE)
  loc <- lwpls_local(fit, test[x_names])
  rw <- loc$training$robust_weight
  expect_gte(sum(order(rw)[1:6] %in% 1:6), 5)
  expect_true(all(loc$training$influence > 0))

  plain <- lwpls_local(lwpls_fit(bad[x_names], bad$y, num_comp = 3), test[x_names])
  expect_true(all(plain$training$robust_weight == 1))
})

test_that("several outcomes and classification", {
  fit <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 2)
  loc <- lwpls_local(fit, test[x_names])
  expect_named(loc$training, c(".row", "y", "y2", "influence", "robust_weight"))
  expect_equal(levels(loc$coefficients$outcome), c("y", "y2"))
  expect_equal(nrow(loc$coefficients), nrow(test) * 6 * 2)
  expect_true(all(c(".pred_y", ".pred_y2") %in% names(loc$reliability)))

  fit_cls <- lwpls_fit(Species ~ ., data = iris, num_comp = 2)
  loc_cls <- lwpls_local(fit_cls, iris[c(1, 51, 101), ])
  expect_equal(as.character(loc_cls$reliability$.pred_class), c("setosa", "versicolor", "virginica"))
  expect_equal(levels(loc_cls$training$.class), levels(iris$Species))
})

test_that("missing values and other interfaces", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2)
  new <- test[1:4, x_names]
  new$x3[2] <- NA
  loc <- lwpls_local(fit, new)
  expect_true(is.na(loc$reliability$status[2]))
  expect_true(all(is.na(loc$coefficients$coefficient[loc$coefficients$.row == 2])))
  expect_equal(nrow(loc$reliability), 4)

  all_na <- new
  all_na$x1 <- NA_real_
  expect_true(all(is.na(lwpls_local(fit, all_na)$reliability$status)))

  spec <- lwpls(num_comp = 2) |> parsnip::set_mode("regression")
  parsnip_fit <- parsnip::fit(spec, y ~ ., data = train[c(x_names, "y")])
  expect_equal(
    lwpls_local(parsnip_fit, test)$reliability,
    lwpls_local(fit, test[x_names])$reliability
  )

  skip_if_not_installed("workflows")
  wflow <- parsnip::fit(workflows::workflow(y ~ ., spec), train[c(x_names, "y")])
  expect_equal(
    lwpls_local(wflow, test)$reliability,
    lwpls_local(fit, test[x_names])$reliability
  )

  expect_snapshot(error = TRUE, lwpls_local(list()))
  expect_snapshot(error = TRUE, lwpls_local(fit, test[x_names], foo = 1))
  expect_snapshot(error = TRUE, lwpls_local(fit, test[x_names], level = 1))
  expect_snapshot(error = TRUE, lwpls_local(fit, test[x_names], limits = "exact"))
  lm_fit <- parsnip::fit(parsnip::linear_reg(), y ~ ., data = train[c(x_names, "y")])
  expect_snapshot(error = TRUE, lwpls_local(lm_fit, test))
})

test_that("spectral axis", {
  expect_equal(spectral_axis(c("x_001", "x_002"))$position, c(1, 2))
  expect_equal(spectral_axis(c("x_001", "x_002"))$label, "Predictor")
  ax <- spectral_axis(c("nm_1650.5", "nm_1652", "nm_1654"))
  expect_equal(ax$position, c(1650.5, 1652, 1654))
  expect_equal(ax$label, "Wavelength")
  expect_equal(spectral_axis(c("a", "b", "c"))$position, 1:3)
  expect_equal(spectral_axis(c("a1", "b1"))$position, 1:2)
  expect_equal(spectral_axis(c("a", "b"), wavelength = c(400, 410))$position, c(400, 410))
  expect_snapshot(error = TRUE, spectral_axis(c("a", "b"), wavelength = 1:3))
})

test_that("print method", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, robust = TRUE, sparsity = 0.3)
  expect_snapshot(lwpls_local(fit, test[x_names]))
})

# Plots -----------------------------------------------------------------------

test_that("plots", {
  fit <- lwpls_fit(train[x_names], train$y, num_comp = 2, robust = TRUE, sparsity = 0.3)
  far <- test[1:2, x_names] + 6
  loc <- lwpls_local(fit, rbind(test[x_names], far), wavelength = seq(400, 900, by = 100))
  builds <- function(plot) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
  }

  rel <- plot_reliability(loc)
  builds(rel)
  expect_equal(rel$labels$x, quote(italic(T)^2 ~ "/" ~ "95% empirical limit"))
  builds(plot_reliability(loc, log = TRUE, labels = 0))
  theory <- plot_reliability(loc, reference = "theoretical")
  builds(theory)
  expect_equal(theory$labels$y, quote(italic(Q) ~ "/" ~ "95% theoretical limit"))
  pct <- plot_reliability(loc, reference = "percentile")
  builds(pct)
  expect_equal(pct$labels$x, expression(italic(T)^2 ~ "percentile rank"))
  expect_snapshot(error = TRUE, plot_reliability(loc, reference = "percentile", log = TRUE))

  coef_plot <- plot_coefficients(loc)
  builds(coef_plot)
  expect_equal(coef_plot$labels$x, "Wavelength")
  builds(plot_coefficients(loc, type = "vip", style = "heatmap", reverse = TRUE))
  builds(plot_coefficients(loc, style = "heatmap", axis_label = "nm"))
  builds(plot_coefficients(loc, reverse = TRUE))

  builds(plot_robust_weights(loc))
  builds(plot_selection(loc))
  builds(plot_selection(loc, by = "model", reverse = TRUE))
  builds(plot_selection(loc, by = "model", spectrum = FALSE))

  for (type in c("reliability", "coefficients", "robust_weights", "selection")) {
    builds(ggplot2::autoplot(loc, type = type))
  }

  few <- lwpls_local(lwpls_fit(train[x_names], train$y, neighbors = 3), test[x_names])
  expect_match(plot_reliability(few)$labels$caption, "queries with too few neighbors")
})

test_that("plots for several outcomes and classes", {
  fit <- lwpls_fit(train[x_names], train[c("y", "y2")], num_comp = 2, robust = TRUE)
  loc <- lwpls_local(fit, test[x_names])
  expect_equal(plot_coefficients(loc, outcome = "y2")$labels$subtitle, "Outcome: y2")
  expect_equal(plot_robust_weights(loc, outcome = "y2")$labels$x, "Reference value (y2)")
  expect_snapshot(error = TRUE, plot_coefficients(loc, outcome = "z"))

  fit_cls <- lwpls_fit(Species ~ ., data = iris, num_comp = 2, robust = TRUE)
  loc_cls <- lwpls_local(fit_cls, iris[c(1, 51, 101), ])
  expect_no_error(ggplot2::ggplot_build(plot_coefficients(loc_cls)))
  expect_no_error(ggplot2::ggplot_build(plot_coefficients(loc_cls, style = "heatmap")))
  expect_no_error(ggplot2::ggplot_build(plot_robust_weights(loc_cls)))

  xy_fit <- parsnip::fit_xy(
    lwpls(num_comp = 2) |> parsnip::set_mode("regression") |> parsnip::set_engine("lwpls", robust = TRUE),
    x = train[x_names], y = train$y
  )
  expect_equal(plot_robust_weights(lwpls_local(xy_fit, test[x_names]))$labels$x, "Reference value")
})

test_that("plots reject unsuitable models", {
  loc <- lwpls_local(lwpls_fit(train[x_names], train$y), test[x_names])
  expect_snapshot(error = TRUE, plot_robust_weights(loc))
  expect_snapshot(error = TRUE, plot_selection(loc))
  expect_snapshot(error = TRUE, plot_reliability(list()))
})
