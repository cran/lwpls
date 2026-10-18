## -----------------------------------------------------------------------------
pkgs <- c("modeldata", "rsample", "tune", "workflows", "ggplot2")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  eval = rlang::is_installed(pkgs)
)

## -----------------------------------------------------------------------------
library(lwpls)
library(parsnip)
library(rsample)
library(tune)
library(workflows)
library(yardstick)
library(ggplot2)

data(meats, package = "modeldata")
meats <- meats[, !names(meats) %in% c("fat", "protein")]

set.seed(123)
meats_split <- initial_split(meats, prop = 0.75)
meats_train <- training(meats_split)
meats_test <- testing(meats_split)

## -----------------------------------------------------------------------------
contaminate <- function(data, fraction, seed = 1) {
  set.seed(seed)
  bad <- sample(nrow(data), round(fraction * nrow(data)))
  shift <- sample(c(-1, 1), length(bad), replace = TRUE) * runif(length(bad), 15, 30)
  data$water[bad] <- data$water[bad] + shift
  data
}

meats_train_20 <- contaminate(meats_train, 0.2)

## -----------------------------------------------------------------------------
set.seed(2)
folds_20 <- vfold_cv(meats_train_20, v = 5)
grid <- expand.grid(num_comp = c(4, 6, 8, 10, 12), localization = c(0.5, 1, 2, 4))

tune_lwpls <- function(robust) {
  wflow <- workflow(
    water ~ .,
    lwpls(num_comp = tune(), localization = tune()) |>
      set_mode("regression") |>
      set_engine("lwpls", robust = robust)
  )
  res <- tune_grid(wflow, folds_20, grid = grid, metrics = metric_set(trimmed_rmse))
  finalize_workflow(wflow, select_best(res, metric = "trimmed_rmse"))
}

lwpls_wflow <- tune_lwpls(robust = FALSE)
robust_wflow <- tune_lwpls(robust = TRUE)

## -----------------------------------------------------------------------------
fractions <- c(0, 0.1, 0.2, 0.3)
test_rmse <- function(wflow, fraction) {
  fitted <- fit(wflow, contaminate(meats_train, fraction))
  rmse_vec(meats_test$water, predict(fitted, meats_test)$.pred)
}

results <- data.frame(
  fraction = rep(fractions, 2),
  model = rep(c("LW-PLS", "Robust LW-PLS"), each = length(fractions)),
  rmse = c(
    vapply(fractions, test_rmse, numeric(1), wflow = lwpls_wflow),
    vapply(fractions, test_rmse, numeric(1), wflow = robust_wflow)
  )
)
results

## -----------------------------------------------------------------------------
ggplot(results, aes(100 * fraction, rmse, color = model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("LW-PLS" = "#969696", "Robust LW-PLS" = "#08519C")) +
  labs(
    x = "Wrong reference values in the training set (%)",
    y = "Test RMSE (water, %)",
    color = NULL
  ) +
  theme_bw()

## -----------------------------------------------------------------------------
set.seed(2)
folds <- vfold_cv(meats_train, v = 5)

sparse_wflow <- workflow(
  water ~ .,
  lwpls(num_comp = tune(), localization = 0.5, sparsity = tune()) |>
    set_mode("regression")
)

sparse_res <- tune_grid(
  sparse_wflow,
  resamples = folds,
  grid = expand.grid(num_comp = 1:15, sparsity = c(0, 0.2, 0.4, 0.6, 0.8)),
  metrics = metric_set(rmse)
)
show_best(sparse_res, metric = "rmse", n = 5)

## -----------------------------------------------------------------------------
autoplot(sparse_res, metric = "rmse") +
  scale_color_viridis_d(option = "mako", end = 0.9) +
  theme_bw()

## -----------------------------------------------------------------------------
sparse_wflow |>
  finalize_workflow(select_best(sparse_res, metric = "rmse")) |>
  last_fit(meats_split) |>
  collect_metrics()

## -----------------------------------------------------------------------------
sprm_fit <-
  lwpls(num_comp = 6, localization = 2, sparsity = 0.4) |>
  set_mode("regression") |>
  set_engine("lwpls", robust = TRUE) |>
  fit(water ~ ., data = meats_train_20)

sprm_fit

rmse_vec(meats_test$water, predict(sprm_fit, meats_test)$.pred)

