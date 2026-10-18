## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.2,
  eval = rlang::is_installed(c("modeldata", "rsample"))
)

## -----------------------------------------------------------------------------
library(lwpls)
library(rsample)

data(meats, package = "modeldata")
meats <- meats[, !names(meats) %in% c("fat", "protein")]

set.seed(123)
meats_split <- initial_split(meats, prop = 0.75)
meats_train <- training(meats_split)
meats_test <- testing(meats_split)

## -----------------------------------------------------------------------------
wavelengths <- seq(850, 1050, length.out = 100)

fit <- lwpls_fit(water ~ ., data = meats_train, num_comp = 10, localization = 0.5)

## -----------------------------------------------------------------------------
drifted <- meats_test[1:3, ]
drifted[, 1:100] <- drifted[, 1:100] + 0.05 * seq(0, 1, length.out = 100)
new_spectra <- rbind(meats_test, drifted)

local <- lwpls_local(fit, new_spectra, wavelength = wavelengths)
local

## -----------------------------------------------------------------------------
plot_reliability(local, log = TRUE)

## -----------------------------------------------------------------------------
rel <- local$reliability
rel$abs_error <- abs(rel$.pred - new_spectra$water)
aggregate(abs_error ~ status, data = rel, FUN = mean)

## -----------------------------------------------------------------------------
theoretical <- lwpls_local(
  fit,
  new_spectra,
  wavelength = wavelengths,
  limits = "theoretical"
)
table(
  empirical = local$reliability$status,
  theoretical = theoretical$reliability$status
)

## -----------------------------------------------------------------------------
plot_reliability(local, reference = "percentile")

## -----------------------------------------------------------------------------
plot_coefficients(local)

## -----------------------------------------------------------------------------
plot_coefficients(local, style = "heatmap")

## -----------------------------------------------------------------------------
set.seed(1)
bad <- sample(nrow(meats_train), round(0.2 * nrow(meats_train)))
contaminated <- meats_train
contaminated$water[bad] <- contaminated$water[bad] +
  sample(c(-1, 1), length(bad), replace = TRUE) * runif(length(bad), 15, 30)

robust_fit <- lwpls_fit(
  water ~ .,
  data = contaminated,
  num_comp = 6,
  localization = 2,
  robust = TRUE
)
robust_local <- lwpls_local(robust_fit, meats_test, wavelength = wavelengths)

## -----------------------------------------------------------------------------
plot_robust_weights(robust_local, labels = 8)

## -----------------------------------------------------------------------------
lowest <- order(robust_local$training$robust_weight)[seq_along(bad)]
sum(lowest %in% bad) / length(bad)

## -----------------------------------------------------------------------------
sparse_fit <- lwpls_fit(
  water ~ .,
  data = meats_train,
  num_comp = 10,
  localization = 0.5,
  sparsity = 0.4
)
sparse_local <- lwpls_local(sparse_fit, meats_test, wavelength = wavelengths)
plot_selection(sparse_local)

## -----------------------------------------------------------------------------
plot_coefficients(local, axis_label = "Wavelength (nm)") +
  ggplot2::labs(title = "Local regression vectors, water content")

## -----------------------------------------------------------------------------
# ggplot2::ggsave("coefficients.pdf", width = 17, height = 10, units = "cm")

