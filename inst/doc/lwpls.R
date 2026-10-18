## -----------------------------------------------------------------------------
pkgs <- c("modeldata", "rsample", "tune", "workflows", "yardstick", "ggplot2")
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

## -----------------------------------------------------------------------------
data(meats, package = "modeldata")
meats <- meats[, !names(meats) %in% c("fat", "protein")]

## -----------------------------------------------------------------------------
set.seed(123)
meats_split <- initial_split(meats, prop = 0.75)
meats_train <- training(meats_split)
meats_test <- testing(meats_split)
meats_folds <- vfold_cv(meats_train, v = 5)

## -----------------------------------------------------------------------------
lwpls_spec <-
  lwpls(num_comp = tune(), localization = tune()) |>
  set_mode("regression")

lwpls_wflow <- workflow(water ~ ., lwpls_spec)

lwpls_grid <- dials::grid_regular(
  dials::num_comp(c(1, 15)),
  localization(),
  levels = c(num_comp = 15, localization = 8)
)

lwpls_res <- tune_grid(
  lwpls_wflow,
  resamples = meats_folds,
  grid = lwpls_grid,
  metrics = metric_set(rmse, rsq)
)

## -----------------------------------------------------------------------------
show_best(lwpls_res, metric = "rmse", n = 5)

## -----------------------------------------------------------------------------
autoplot(lwpls_res, metric = "rmse") +
  scale_color_viridis_d(option = "mako", end = 0.9) +
  theme_bw()

## -----------------------------------------------------------------------------
# rmse_grid <- collect_metrics(lwpls_res)
# rmse_grid <- rmse_grid[rmse_grid$.metric == "rmse", ]
# best <- rmse_grid[which.min(rmse_grid$mean), ]
# locs <- sort(unique(rmse_grid$localization))
# # rows: localization, columns: num_comp
# rmse_mat <- with(rmse_grid, tapply(mean, list(localization, num_comp), identity))
# hover <- "num_comp: %{x}<br>localization: %{y:.4g}<br>RMSE: %{z:.3f}<extra></extra>"
# 
# plotly::plot_ly() |>
#   plotly::add_surface(
#     x = sort(unique(rmse_grid$num_comp)),
#     y = locs,
#     z = rmse_mat,
#     colorscale = list(list(0, "#cde2fb"), list(0.5, "#3987e5"), list(1, "#0d366b")),
#     showscale = FALSE, hovertemplate = hover
#     ) |>
#   plotly::add_markers(
#     x = best$num_comp,
#     y = best$localization,
#     z = best$mean,
#     marker = list(size = 6, color = "#eb6834", line = list(color = "white", width = 2)),
#     hovertemplate = paste0("<b>Best</b><br>", hover), showlegend = FALSE
#     ) |>
#   plotly::layout(
#     scene = list(
#       xaxis = list(title = "num_comp"),
#       yaxis = list(
#         title = "localization", type = "log",
#         tickvals = locs, ticktext = paste0("2<sup>", log2(locs), "</sup>")),
#       zaxis = list(title = "RMSE"),
#       camera = list(eye = list(x = 1.7, y = 1.6, z = 0.8))),
#     margin = list(l = 0, r = 0, t = 0, b = 0)
#   )

## -----------------------------------------------------------------------------
best_params <- select_best(lwpls_res, metric = "rmse")
best_params

lwpls_final <-
  lwpls_wflow |>
  finalize_workflow(best_params) |>
  last_fit(meats_split)

collect_metrics(lwpls_final)

## -----------------------------------------------------------------------------
pls_wflow <- workflow(
  water ~ .,
  lwpls(num_comp = tune(), localization = 1e6) |> set_mode("regression")
)
pls_res <- tune_grid(
  pls_wflow,
  resamples = meats_folds,
  grid = data.frame(num_comp = 1:15),
  metrics = metric_set(rmse, rsq)
)

pls_final <-
  pls_wflow |>
  finalize_workflow(select_best(pls_res, metric = "rmse")) |>
  last_fit(meats_split)

collect_metrics(pls_final)

## -----------------------------------------------------------------------------
test_preds <- rbind(
  data.frame(model = "Global PLS", collect_predictions(pls_final)),
  data.frame(model = "LW-PLS", collect_predictions(lwpls_final))
)

ggplot(test_preds, aes(water, .pred)) +
  geom_abline(color = "grey50", linetype = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~model) +
  coord_obs_pred() +
  labs(x = "Observed water (%)", y = "Predicted water (%)") +
  theme_bw()

## -----------------------------------------------------------------------------
similarity_wflow <- workflow(
  water ~ .,
  lwpls(num_comp = tune(), localization = tune(), similarity = tune()) |>
    set_mode("regression")
)

similarity_grid <- dials::grid_regular(
  dials::num_comp(c(1, 15)),
  localization(),
  similarity(),
  levels = c(num_comp = 15, localization = 8, similarity = 2)
)

similarity_res <- tune_grid(
  similarity_wflow,
  resamples = meats_folds,
  grid = similarity_grid,
  metrics = metric_set(rmse)
)

# Best configuration of each similarity index
similarity_perf <- collect_metrics(similarity_res)
best_by_index <- lapply(
  split(similarity_perf, similarity_perf$similarity),
  function(x) x[which.min(x$mean), ]
)
do.call(rbind, best_by_index)[, c("similarity", "num_comp", "localization", "mean", "std_err")]

## -----------------------------------------------------------------------------
data(parabolic, package = "modeldata")

set.seed(91)
para_split <- initial_split(parabolic, strata = class)
para_train <- training(para_split)
para_test <- testing(para_split)
para_folds <- vfold_cv(para_train, v = 5, strata = class)

## -----------------------------------------------------------------------------
lwplsda_wflow <- workflow(
  class ~ .,
  lwpls(
    num_comp = 2, 
    localization = tune()
    ) |>
    set_mode("classification")
  )

lwplsda_res <- tune_grid(
  lwplsda_wflow,
  resamples = para_folds,
  grid = data.frame(localization = 2^seq(-6, 5)),
  metrics = metric_set(roc_auc, accuracy)
)

show_best(lwplsda_res, metric = "roc_auc", n = 3)

## -----------------------------------------------------------------------------
autoplot(lwplsda_res) +
  theme_bw()

## -----------------------------------------------------------------------------
lwplsda_final <-
  lwplsda_wflow |>
  finalize_workflow(select_best(lwplsda_res, metric = "roc_auc")) |>
  last_fit(para_split, metrics = metric_set(roc_auc, accuracy))

plsda_final <-
  workflow(
    class ~ .,
    lwpls(num_comp = 2, localization = 1e6) |> set_mode("classification")) |>
  last_fit(para_split, metrics = metric_set(roc_auc, accuracy))

rbind(
  data.frame(model = "LW-PLS-DA", collect_metrics(lwplsda_final)),
  data.frame(model = "Global PLS-DA", collect_metrics(plsda_final))
)[, c("model", ".metric", ".estimate")]

## -----------------------------------------------------------------------------
para_grid <- expand.grid(
  X1 = seq(min(parabolic$X1), max(parabolic$X1), length.out = 100),
  X2 = seq(min(parabolic$X2), max(parabolic$X2), length.out = 100)
)

region_probs <- function(final, model) {
  probs <- predict(extract_workflow(final), para_grid, type = "prob")
  data.frame(para_grid, model = model, .pred_Class1 = probs$.pred_Class1)
}

regions <- rbind(
  region_probs(plsda_final, "Global PLS-DA"),
  region_probs(lwplsda_final, "LW-PLS-DA")
)

ggplot(regions, aes(X1, X2)) +
  geom_raster(aes(fill = .pred_Class1), alpha = 0.8) +
  geom_contour(aes(z = .pred_Class1), breaks = 0.5, color = "black") +
  geom_point(data = para_test, aes(shape = class), color = "white", size = 2.6) +
  geom_point(data = para_test, aes(shape = class, color = class), size = 1.6) +
  facet_wrap(~model) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1)) +
  scale_color_manual(values = c(Class1 = "#08306B", Class2 = "#67000D")) +
  coord_equal(expand = FALSE) +
  labs(fill = "Probability\nof Class1", color = NULL, shape = NULL) +
  theme_bw()

## -----------------------------------------------------------------------------
test_roc <- function(final, model) {
  roc <- roc_curve(collect_predictions(final), truth = class, .pred_Class1)
  data.frame(model = model, roc)
}

roc_data <- rbind(
  test_roc(lwplsda_final, "LW-PLS-DA"),
  test_roc(plsda_final, "Global PLS-DA")
)

ggplot(roc_data, aes(1 - specificity, sensitivity, color = model)) +
  geom_abline(color = "grey50", linetype = 2) +
  geom_path(linewidth = 0.8) +
  scale_color_manual(values = c("LW-PLS-DA" = "#08519C", "Global PLS-DA" = "#969696")) +
  coord_equal() +
  labs(color = NULL) +
  theme_bw()

## -----------------------------------------------------------------------------
collect_predictions(lwplsda_final) |>
  conf_mat(truth = class, estimate = .pred_class) |>
  autoplot(type = "heatmap")

## -----------------------------------------------------------------------------
fit <- lwpls_fit(water ~ ., data = meats_train, num_comp = 10, localization = 0.5)
fit

predict(fit, head(meats_test))
predict(fit, head(meats_test), num_comp = 3)

## -----------------------------------------------------------------------------
parsnip_fit <- extract_fit_parsnip(lwpls_final)
multi_predict(parsnip_fit, head(meats_test), num_comp = c(5, 10, 15))$.pred[[1]]

## -----------------------------------------------------------------------------
lwpls(num_comp = 10, localization = 0.5) |>
  set_mode("regression") |>
  set_engine("lwpls", scale = FALSE) |>
  translate()

