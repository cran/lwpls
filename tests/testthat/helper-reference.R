# Direct R transcription of `lwpls()` from https://github.com/hkaneko1985/lwpls
# (Copyright (c) 2018 Hiromasa Kaneko, MIT License) (NIPALS with explicit
# deflation). `x_train`, `y_train` and `x_test` must be
# autoscaled with the training statistics. Returns an
# nrow(x_test) x max_component_number matrix of predictions.
#
# `gamma` (optional) is the p x r projection used to measure distances,
# d_n = ||gamma' (x_n - x_q)||: Gamma = X'y / ||X'y|| gives covariance-based
# LW-PLS (Hazama & Kano 2015, eqs. 28 and 33).
kaneko_lwpls <- function(x_train, y_train, x_test, max_component_number,
                         lambda_in_similarity, gamma = NULL) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- matrix(y_train, ncol = 1)
  if (is.null(gamma)) {
    gamma <- diag(ncol(x_train))
  }
  est <- matrix(0, nrow(x_test), max_component_number)
  for (i in seq_len(nrow(x_test))) {
    xq <- x_test[i, , drop = FALSE]
    d <- sqrt(colSums((crossprod(gamma, t(x_train) - drop(xq)))^2))
    s <- exp(-d / stats::sd(d) / lambda_in_similarity)
    yw <- sum(s * y_train) / sum(s)
    xw <- colSums(s * x_train) / sum(s)
    cy <- y_train - yw
    cx <- sweep(x_train, 2, xw)
    cq <- xq - xw
    est[i, ] <- yw
    for (a in seq_len(max_component_number)) {
      wa <- crossprod(cx, s * cy)
      wa <- wa / sqrt(sum(wa^2))
      ta <- cx %*% wa
      pa <- crossprod(cx, s * ta) / sum(s * ta^2)
      qa <- sum(s * cy * ta) / sum(s * ta^2)
      tqa <- drop(cq %*% wa)
      est[i, a:max_component_number] <- est[i, a:max_component_number] + tqa * qa
      cx <- cx - ta %*% t(pa)
      cy <- cy - ta * qa
      cq <- cq - tqa * t(pa)
    }
  }
  est
}

# Kaneko's LW-PLS on the original scale (autoscaling as in his demo script).
kaneko_predict <- function(x_train, y_train, x_test, num_comp, localization,
                           scale = TRUE, covariance = FALSE) {
  mu <- colMeans(x_train)
  s <- if (scale) apply(x_train, 2, stats::sd) else rep(1, ncol(x_train))
  y_mu <- mean(y_train)
  y_s <- if (scale) stats::sd(y_train) else 1
  xs <- scale(x_train, mu, s)
  ys <- (y_train - y_mu) / y_s
  gamma <- NULL
  if (covariance) {
    gamma <- crossprod(xs, ys)
    gamma <- gamma / sqrt(sum(gamma^2))
  }
  est <- kaneko_lwpls(
    xs,
    ys,
    scale(x_test, mu, s),
    num_comp,
    localization,
    gamma = gamma
  )
  est * y_s + y_mu
}

# Simulated nonlinear regression data.
sim_reg <- function(n = 120, p = 6, seed = 1) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
    y <- sin(x[, 1]) + x[, 2]^2 + 0.5 * x[, 3] + stats::rnorm(n, sd = 0.1)
    y2 <- x[, 4] * x[, 5] + stats::rnorm(n, sd = 0.1)
  })
  dat <- as.data.frame(x)
  dat$y <- y
  dat$y2 <- y2
  dat
}

# ------------------------------------------------------------------------------
# Sparse and robust local models
#
# Plain R implementations of weighted SNIPLS (Hoffmann et al. 2015) and local
# partial robust M-regression (Serneels et al. 2005). Checked outside the test
# suite against sprm::snipls() (identical except where sprm keeps sign-flipped
# weights for previously selected variables) and chemometrics::prm()
# (identical, to 1e-14, once its scores are centered per column: prm() centers
# them with `T - mt`, which recycles the center down the rows when a > 1).

# Weighted median; equals stats::median() for equal weights.
wmedian <- function(x, w) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w)
  half <- sum(w) / 2
  tol <- 1e-9 * sum(w)
  k <- which(cw >= half - tol)[1]
  if (abs(cw[k] - half) <= tol && k < length(x)) (x[k] + x[k + 1]) / 2 else x[k]
}

# Weighted L1-median (spatial median) by Weiszfeld iterations.
wl1median <- function(X, w, tol = 1e-12, maxit = 1000) {
  if (ncol(X) == 1) return(wmedian(X[, 1], w))
  m <- apply(X, 2, wmedian, w = w)
  for (it in seq_len(maxit)) {
    d <- sqrt(rowSums(sweep(X, 2, m)^2))
    d <- pmax(d, 1e-12 * max(d, 1e-300))
    v <- w / d
    m_new <- colSums(v * X) / sum(v)
    step <- sqrt(sum((m_new - m)^2))
    m <- m_new
    if (step <= tol * max(1, sqrt(sum(m^2)))) break
  }
  m
}

fair_weight <- function(z, c) 1 / (1 + abs(z / c))^2

hampel_weight <- function(z, k) {
  z <- abs(z)
  w <- rep(1, length(z))
  i <- z > k[1] & z <= k[2]
  w[i] <- k[1] / z[i]
  i <- z > k[2] & z <= k[3]
  w[i] <- k[1] * (k[3] - z[i]) / ((k[3] - k[2]) * z[i])
  w[z > k[3]] <- 0
  w
}

# Weighted NIPALS on centered data with optional SNIPLS sparsity. Returns the
# training scores, Y loadings, t'Wt and the query predictions (centered) for
# 1..a components.
wnipals <- function(Xc, Yc, w, a, sparsity = 0, xq = NULL) {
  p <- ncol(Xc)
  q <- ncol(Yc)
  Xa <- Xc
  Ya <- Yc
  xqa <- xq
  active <- rep(FALSE, p)
  Tm <- matrix(0, nrow(Xc), a)
  C <- matrix(0, q, a)
  tt <- numeric(a)
  yq <- matrix(0, q, a)
  cur <- rep(0, q)
  for (j in seq_len(a)) {
    M <- crossprod(Xa, w * Ya)
    wv <- if (q == 1) M[, 1] else svd(M, nu = 1, nv = 0)$u[, 1]
    wv <- wv / sqrt(sum(wv^2))
    if (sparsity > 0) {
      thr <- sparsity * max(abs(wv))
      active <- active | (abs(wv) >= thr)
      wv <- sign(wv) * pmax(abs(wv) - thr, 0)
      wv <- wv / sqrt(sum(wv^2))
    }
    t <- drop(Xa %*% wv)
    ttj <- sum(w * t^2)
    pv <- drop(crossprod(Xa, w * t)) / ttj
    if (sparsity > 0) pv[!active] <- 0
    cv <- drop(crossprod(Ya, w * t)) / ttj
    if (!is.null(xqa)) {
      tq <- sum(xqa * wv)
      cur <- cur + tq * cv
      xqa <- xqa - tq * pv
    }
    yq[, j] <- cur
    Xa <- Xa - tcrossprod(t, pv)
    Ya <- Ya - tcrossprod(t, cv)
    Tm[, j] <- t
    C[, j] <- cv
    tt[j] <- ttj
  }
  list(T = Tm, C = C, tt = tt, yq = yq)
}

# Standardized distances of the rows of a residual matrix and their weights.
residual_weights <- function(R, om, fun, fair_c, probs, df) {
  R <- as.matrix(R)
  rc <- sweep(R, 2, apply(R, 2, wmedian, w = om))
  s <- apply(abs(rc), 2, wmedian, w = om)
  if (fun == "fair") {
    d <- sqrt(rowSums(sweep(rc, 2, s, "/")^2))
    fair_weight(d / wmedian(d, om), fair_c)
  } else {
    d <- sqrt(rowSums(sweep(rc, 2, 1.4826 * s, "/")^2))
    hampel_weight(d, sqrt(qchisq(probs, df)))
  }
}

leverage_weights <- function(Tm, tt, om, fun, fair_c, probs) {
  if (fun == "fair") {
    Tn <- sweep(Tm, 2, sqrt(tt), "/")
    d <- sqrt(rowSums(sweep(Tn, 2, wl1median(Tn, om))^2))
    fair_weight(d / wmedian(d, om), fair_c)
  } else {
    tc <- sweep(Tm, 2, apply(Tm, 2, wmedian, w = om))
    s <- 1.4826 * apply(abs(tc), 2, wmedian, w = om)
    d <- sqrt(rowSums(sweep(tc, 2, s, "/")^2))
    hampel_weight(d, sqrt(qchisq(probs, ncol(Tm))))
  }
}

# Local PRM prediction of one query with a components.
local_prm <- function(X, Y, xq, omega, a, sparsity = 0, fun = "fair", fair_c = 4,
                      probs = c(0.95, 0.975, 0.999), max_iter = 30, tol = 0.01,
                      classification = FALSE) {
  Y <- as.matrix(Y)
  keep <- omega > 0
  X <- X[keep, , drop = FALSE]
  Y <- Y[keep, , drop = FALSE]
  om <- omega[keep]
  q <- ncol(Y)
  cx <- wl1median(X, om)
  cy <- if (classification) colSums(om * Y) / sum(om) else apply(Y, 2, wmedian, w = om)
  Xc <- sweep(X, 2, cx)
  Yc <- sweep(Y, 2, cy)
  xqc <- xq - cx
  df_res <- if (classification) q - 1 else q

  dx <- sqrt(rowSums(Xc^2))
  wx <- if (fun == "fair") fair_weight(dx / wmedian(dx, om), fair_c) else
    hampel_weight(dx / wmedian(dx, om), qnorm(probs))
  wy <- if (classification) rep(1, length(om)) else
    residual_weights(Yc, om, fun, fair_c, probs, df_res)
  wr <- wx * wy
  wr[wr < 1e-6] <- 1e-6

  gamma <- 1e5
  diff <- 1
  iter <- 1
  while (diff > tol && iter <= max_iter) {
    wt <- om * wr
    fit <- wnipals(Xc, Yc, wt, a, sparsity, xqc)
    gamma_old <- gamma
    gamma <- sqrt(sum(sweep(fit$C, 2, sqrt(fit$tt), "*")^2))
    diff <- abs(gamma - gamma_old) / gamma
    E <- Yc - fit$T %*% t(fit$C)
    wres <- residual_weights(E, om, fun, fair_c, probs, df_res)
    wlev <- leverage_weights(fit$T, fit$tt, om, fun, fair_c, probs)
    wr <- wres * wlev
    wr[wr < 1e-6] <- 1e-6
    iter <- iter + 1
  }
  adj <- if (classification) colSums(wt * E) / sum(wt) else apply(E, 2, wmedian, w = om)
  cy + adj + fit$yq[, a]
}

# Similarity weights as computed by the C++ kernel.
ref_omega <- function(d, localization, robust) {
  s <- stats::sd(d)
  if (robust) {
    mad_d <- 1.4826 * stats::median(abs(d - stats::median(d)))
    if (mad_d > 0) s <- mad_d
  }
  w <- exp(-(d - min(d)) / (s * localization))
  w[w < .Machine$double.eps] <- 0
  w
}

# Predictions of a fitted sparse and/or robust model by the reference code.
ref_predict <- function(fit, new_x, num_comp) {
  qn <- scale(as.matrix(new_x), fit$x_center, fit$x_scale)
  fun <- fit$weight_function
  out <- t(apply(qn, 1, function(xq) {
    d <- sqrt(colSums((t(fit$x) - xq)^2))
    om <- ref_omega(d, fit$localization, fit$robust)
    if (!fit$robust) {
      keep <- om > 0
      X <- fit$x[keep, , drop = FALSE]
      Y <- fit$y[keep, , drop = FALSE]
      o <- om[keep]
      cx <- colSums(o * X) / sum(o)
      cy <- colSums(o * Y) / sum(o)
      pred <- cy + wnipals(sweep(X, 2, cx), sweep(Y, 2, cy), o, num_comp,
        fit$sparsity, xq - cx)$yq[, num_comp]
    } else {
      pred <- local_prm(
        fit$x, fit$y, xq, om, num_comp,
        sparsity = fit$sparsity,
        fun = fun,
        fair_c = if (fun == "fair") fit$robust_constant else 4,
        probs = if (fun == "hampel") fit$robust_constant else c(0.95, 0.975, 0.999),
        max_iter = fit$max_iter,
        classification = fit$mode == "classification"
      )
    }
    pred * fit$y_scale + fit$y_center
  }))
  out <- unname(out)
  if (ncol(fit$y) == 1) drop(out) else out
}
