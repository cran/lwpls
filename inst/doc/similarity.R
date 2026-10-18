## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3.2,
  fig.retina = 2,
  fig.align = "center",
  out.width = "100%"
)

## -----------------------------------------------------------------------------
library(lwpls)
library(ggplot2)

## -----------------------------------------------------------------------------
theme_article <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      panel.border = element_rect(colour = "grey35", linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", colour = "grey35",
                                       linewidth = 0.4),
      axis.ticks = element_line(colour = "grey35", linewidth = 0.3),
      axis.text = element_text(colour = "grey20"),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

index_colors <- c("Euclidean" = "#2a78d6", "Covariance-based" = "#eb6834")
weight_colors <- c("white", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab",
                   "#0d366b")

index_label <- function(similarity) {
  factor(ifelse(similarity == "euclidean", "Euclidean", "Covariance-based"),
         levels = names(index_colors))
}

# Test RMSE versus localization, with the best value of each index and the
# global PLS model for reference.
plot_rmse <- function(curves, global) {
  curves$index <- index_label(curves$similarity)
  best <- do.call(rbind, lapply(split(curves, curves$index), function(d) {
    d[which.min(d$rmse), ]
  }))
  best$label <- sprintf("%.3f", best$rmse)
  ggplot(curves, aes(localization, rmse, colour = index)) +
    geom_hline(yintercept = global, linetype = "dashed", colour = "grey35",
               linewidth = 0.4) +
    annotate("text", x = 1, y = global, vjust = -0.5, size = 3,
             colour = "grey25", label = sprintf("global PLS (%.3f)", global)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.1) +
    geom_point(data = best, size = 3, shape = 21, fill = "white",
               stroke = 0.9) +
    geom_text(data = best, aes(label = label), vjust = 2, size = 3,
              show.legend = FALSE) +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = c(0.14, 0.08))) +
    scale_colour_manual(values = index_colors, name = NULL) +
    labs(x = expression("Localization " * varphi), y = "Test RMSE") +
    theme_article() +
    theme(legend.position = "inside", legend.position.inside = c(0.985, 0.03),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = NA))
}

## -----------------------------------------------------------------------------
kernel <- expand.grid(u = seq(0, 4, length.out = 400), phi = c(0.25, 0.5, 1, 2))
kernel$weight <- exp(-kernel$u / kernel$phi)
crossings <- data.frame(u = c(0.25, 0.5, 1, 2), phi = c(0.25, 0.5, 1, 2),
                        weight = exp(-1))

ggplot(kernel, aes(u, weight, colour = factor(phi))) +
  geom_hline(yintercept = exp(-1), linetype = "dotted", colour = "grey45",
             linewidth = 0.35) +
  geom_line(linewidth = 0.7) +
  geom_point(data = crossings, size = 1.8) +
  scale_colour_manual(values = weight_colors[3:6], name = expression(varphi)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = expression(d[i] / sigma[d]), y = expression(omega[i])) +
  theme_article() +
  theme(legend.position = "inside", legend.position.inside = c(0.88, 0.68))

## -----------------------------------------------------------------------------
set.seed(2026)
n <- 250
x1 <- rnorm(n)
x2 <- 0.75 * x1 + sqrt(1 - 0.75^2) * rnorm(n)
y <- sin(1.5 * x1) + 0.3 * x2 + rnorm(n, sd = 0.1)
x <- scale(cbind(x1, x2))

gamma <- crossprod(x, y - mean(y))
gamma <- drop(gamma / sqrt(sum(gamma^2)))
gamma

query <- c(0.6, 0.2)
distance <- list(
  "Euclidean" = function(z) sqrt(rowSums(sweep(z, 2, query)^2)),
  "Mahalanobis" = function(z) sqrt(mahalanobis(z, query, cov(x))),
  "Covariance-based" = function(z) abs(drop(sweep(z, 2, query) %*% gamma))
)

## -----------------------------------------------------------------------------
phi <- 0.5
grid_values <- seq(-3, 3, length.out = 241)
grid <- as.matrix(expand.grid(x1 = grid_values, x2 = grid_values))
panels <- names(distance)

fields <- do.call(rbind, lapply(panels, function(nm) {
  d_train <- distance[[nm]](x)
  d_grid <- distance[[nm]](grid)
  data.frame(grid, index = nm,
             weight = pmin(1, exp(-(d_grid - min(d_train)) / (sd(d_train) * phi))))
}))
# Facet labels with the effective number of neighbors of each distance.
facet_labels <- vapply(panels, function(nm) {
  d_train <- distance[[nm]](x)
  w <- exp(-(d_train - min(d_train)) / (sd(d_train) * phi))
  sprintf('"%s" ~~ (italic(n)[eff] == %.0f)', nm, sum(w)^2 / sum(w^2))
}, character(1))
fields$index <- factor(facet_labels[fields$index], levels = facet_labels)
arrow_df <- data.frame(index = factor(facet_labels[["Covariance-based"]],
                                      levels = facet_labels),
                       x = query[1], y = query[2],
                       xend = query[1] + 1.5 * gamma[1],
                       yend = query[2] + 1.5 * gamma[2])
training <- data.frame(x1 = x[, 1], x2 = x[, 2])

ggplot(fields, aes(x1, x2)) +
  geom_raster(aes(fill = weight), interpolate = TRUE) +
  geom_contour(aes(z = weight), breaks = c(0.1, 0.5), colour = "white",
               linewidth = 0.3) +
  geom_point(data = training, shape = 21, size = 0.9, fill = "grey15",
             colour = "white", stroke = 0.15) +
  geom_segment(data = arrow_df, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.6, colour = "black",
               arrow = grid::arrow(length = grid::unit(1.8, "mm"), type = "closed")) +
  geom_text(data = arrow_df, aes(x = xend, y = yend), label = "Gamma",
            parse = TRUE, hjust = -0.3, vjust = -0.1, size = 3.6) +
  annotate("point", x = query[1], y = query[2], shape = 21, size = 2.6,
           fill = "white", colour = "black", stroke = 0.7) +
  facet_wrap(~index, nrow = 1, labeller = label_parsed) +
  scale_fill_gradientn(colours = weight_colors, limits = c(0, 1),
                       name = expression(omega[i] / max(omega))) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  coord_equal(xlim = c(-3, 3), ylim = c(-3, 3), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  labs(x = expression(x[1] ~ "(standardized)"), y = expression(x[2] ~ "(standardized)")) +
  theme_article() +
  theme(legend.position = "bottom", legend.key.width = grid::unit(1.6, "cm"),
        legend.key.height = grid::unit(0.25, "cm"), panel.grid = element_blank(),
        panel.spacing = grid::unit(0.8, "lines"),
        legend.margin = margin(0, 0, 0, 0))

## -----------------------------------------------------------------------------
wavelength <- seq(1100, 1700, by = 5)
band <- function(center, width) exp(-0.5 * ((wavelength - center) / width)^2)

simulate_spectra <- function(n) {
  concentration <- runif(n, 0.05, 1)
  absorbance <- 1 - exp(-2.5 * concentration)
  interferent <- runif(n, 0, 1.5)
  offset <- rnorm(n, sd = 0.15)
  slope <- rnorm(n, sd = 0.1)
  x <- outer(0.8 * absorbance, band(1460, 22)) +
    outer(interferent, band(1330, 70)) +
    outer(offset, rep(1, length(wavelength))) +
    outer(slope, (wavelength - 1400) / 300) +
    matrix(rnorm(n * length(wavelength), sd = 0.004), n)
  colnames(x) <- paste0("nm", wavelength)
  list(x = x, y = concentration)
}

set.seed(7)
train <- simulate_spectra(300)
test <- simulate_spectra(300)

## -----------------------------------------------------------------------------
x_train <- scale(train$x)
gamma_spectra <- crossprod(x_train, train$y - mean(train$y))
gamma_spectra <- drop(gamma_spectra / sqrt(sum(gamma_spectra^2)))
cor(drop(x_train %*% gamma_spectra), train$y)

## -----------------------------------------------------------------------------
shown <- 1:40
spectra_df <- data.frame(
  wavelength = rep(wavelength, each = length(shown)),
  value = as.vector(train$x[shown, ]),
  y = rep(train$y[shown], times = length(wavelength)),
  id = rep(shown, times = length(wavelength)),
  panel = "Absorbance"
)
gamma_df <- data.frame(wavelength = wavelength, value = gamma_spectra,
                       panel = "Gamma[j]")
band_df <- data.frame(wavelength = wavelength,
                      value = band(1460, 22) * max(gamma_spectra),
                      panel = "Gamma[j]")
labels_df <- data.frame(wavelength = c(1330, 1460),
                        value = max(spectra_df$value) * 1.07,
                        label = c("interferent", "analyte"),
                        panel = "Absorbance")
as_panel <- function(df) {
  df$panel <- factor(df$panel, levels = c("Absorbance", "Gamma[j]"))
  df
}
spectra_df <- as_panel(spectra_df)
gamma_df <- as_panel(gamma_df)
band_df <- as_panel(band_df)
labels_df <- as_panel(labels_df)

ggplot(mapping = aes(wavelength, value)) +
  geom_area(data = band_df, fill = "grey88") +
  geom_hline(data = gamma_df[1, ], aes(yintercept = 0), colour = "grey35",
             linewidth = 0.3) +
  geom_line(data = spectra_df, aes(group = id, colour = y), linewidth = 0.35) +
  geom_line(data = gamma_df, colour = index_colors[["Covariance-based"]],
            linewidth = 0.7) +
  geom_text(data = labels_df, aes(label = label), size = 3, colour = "grey25") +
  facet_grid(panel ~ ., scales = "free_y", switch = "y",
             labeller = label_parsed) +
  scale_colour_gradientn(colours = weight_colors[-1], name = "Concentration") +
  labs(x = "Wavelength (nm)", y = NULL) +
  theme_article() +
  theme(strip.placement = "outside", strip.background.y = element_blank(),
        strip.text.y.left = element_text(angle = 90, size = rel(1.1)),
        legend.title = element_text(size = rel(0.9)),
        legend.key.height = grid::unit(0.9, "cm"),
        legend.key.width = grid::unit(0.25, "cm"))

## -----------------------------------------------------------------------------
rmse_by_localization <- function(train, test, num_comp, localization) {
  grid <- expand.grid(localization = localization,
                      similarity = c("euclidean", "covariance"),
                      stringsAsFactors = FALSE)
  grid$rmse <- mapply(function(phi, index) {
    fit <- lwpls_fit(train$x, train$y, num_comp = num_comp,
                     localization = phi, similarity = index)
    sqrt(mean((test$y - predict(fit, test$x)$.pred)^2))
  }, grid$localization, grid$similarity)
  grid
}

global_rmse <- function(train, test, num_comp) {
  fit <- lwpls_fit(train$x, train$y, num_comp = num_comp, localization = 1e6)
  sqrt(mean((test$y - predict(fit, test$x)$.pred)^2))
}

localization_grid <- 10^seq(-1.5, 1, length.out = 14)
spectra_rmse <- rmse_by_localization(train, test, num_comp = 4, localization_grid)
spectra_global <- global_rmse(train, test, num_comp = 4)

# Best localization of each index
do.call(rbind, lapply(split(spectra_rmse, spectra_rmse$similarity),
                      function(d) d[which.min(d$rmse), ]))
spectra_global

## -----------------------------------------------------------------------------
plot_rmse(spectra_rmse, spectra_global)

## -----------------------------------------------------------------------------
query_id <- which.min(abs(test$y - 0.5))
query_x <- test$x[query_id, , drop = FALSE]
z_query <- (drop(query_x) - attr(x_train, "scaled:center")) /
  attr(x_train, "scaled:scale")

d_euclidean <- sqrt(rowSums(sweep(x_train, 2, z_query)^2))
d_covariance <- abs(drop(sweep(x_train, 2, z_query) %*% gamma_spectra))
relative_weights <- function(d, phi) exp(-(d - min(d)) / (sd(d) * phi))

query_fits <- lapply(c(euclidean = "euclidean", covariance = "covariance"),
                     function(index) {
  lwpls_fit(train$x, train$y, num_comp = 4, localization = 0.2,
            similarity = index)
})
sapply(query_fits, function(fit) lwpls_local(fit, query_x)$reliability$n_eff)

## -----------------------------------------------------------------------------
weights_df <- rbind(
  data.frame(similarity = "euclidean", y = train$y,
             weight = relative_weights(d_euclidean, 0.2)),
  data.frame(similarity = "covariance", y = train$y,
             weight = relative_weights(d_covariance, 0.2))
)
weights_df$index <- index_label(weights_df$similarity)

ggplot(weights_df, aes(y, weight, colour = index)) +
  geom_vline(xintercept = test$y[query_id], linetype = "dashed",
             colour = "grey35", linewidth = 0.4) +
  geom_point(size = 0.9, alpha = 0.75) +
  facet_wrap(~index, nrow = 1) +
  scale_colour_manual(values = index_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Concentration of the training sample",
       y = expression(omega[i] / max(omega))) +
  theme_article()

## -----------------------------------------------------------------------------
simulate_symmetric <- function(n) {
  x <- matrix(rnorm(n * 5), n, dimnames = list(NULL, paste0("x", 1:5)))
  list(x = x, y = x[, 1]^2 + rnorm(n, sd = 0.2))
}

set.seed(11)
train_sym <- simulate_symmetric(300)
test_sym <- simulate_symmetric(300)

x_sym <- scale(train_sym$x)
gamma_sym <- crossprod(x_sym, train_sym$y - mean(train_sym$y))
gamma_sym <- drop(gamma_sym / sqrt(sum(gamma_sym^2)))
round(gamma_sym, 2)
cor(drop(x_sym %*% gamma_sym), train_sym$y)

## -----------------------------------------------------------------------------
symmetric_df <- rbind(
  data.frame(panel = "x[1]", coordinate = x_sym[, 1], y = train_sym$y),
  data.frame(panel = "t == X * Gamma", coordinate = drop(x_sym %*% gamma_sym),
             y = train_sym$y)
)
symmetric_df$panel <- factor(symmetric_df$panel,
                             levels = c("x[1]", "t == X * Gamma"))

ggplot(symmetric_df, aes(coordinate, y)) +
  geom_point(size = 0.8, alpha = 0.6, colour = "grey20") +
  facet_wrap(~panel, nrow = 1, labeller = label_parsed) +
  labs(x = "Standardized coordinate", y = "y") +
  theme_article()

## -----------------------------------------------------------------------------
symmetric_rmse <- rmse_by_localization(train_sym, test_sym, num_comp = 2,
                                       localization_grid)
symmetric_global <- global_rmse(train_sym, test_sym, num_comp = 2)
do.call(rbind, lapply(split(symmetric_rmse, symmetric_rmse$similarity),
                      function(d) d[which.min(d$rmse), ]))
symmetric_global

## -----------------------------------------------------------------------------
plot_rmse(symmetric_rmse, symmetric_global)

