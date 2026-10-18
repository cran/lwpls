// Locally-weighted partial least squares (LW-PLS): prediction kernel.
//
// For every query sample x_q, each training sample i receives the similarity
// weight (Kim et al. 2011; Kaneko)
//
//     w_i = exp(-d_i / (sd(d) * localization)),   d_i = ||z_i - z_q||,
//
// where z are the coordinates used to measure similarity: the predictors
// themselves (Euclidean distance) or their projection on the covariance
// direction X'Y (CbLW-PLS, Hazama & Kano 2015). A weighted PLS model is then
// fitted around the query. The weighted PLS uses
// the kernel formulation of Dayal & MacGregor (1997): only the p x q
// cross-product matrix X'WY is deflated and the scores are computed through
// the R = W(P'W)^-1 weights, so the n x p training matrix is never copied or
// deflated. Queries are processed in blocks so that the expensive steps
// (distances, local means, scores and loadings) are matrix-matrix products
// that run through level-3 BLAS.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace {

// Weights smaller than this (relative to the nearest sample, whose weight is
// one) are set to zero. They have no measurable influence on the local model
// but would otherwise let numerical noise dominate when the model is very
// local.
const double weight_floor = std::numeric_limits<double>::epsilon();

// A local covariance ||X_c' W Y_c|| below this fraction of its Cauchy-Schwarz
// bound sqrt(SS_x SS_y) is rounding noise: no component is extracted.
const double cov_floor = std::sqrt(std::numeric_limits<double>::epsilon());

// Median of the values in `v` (reordered in place).
double median_inplace(std::vector<double>& v) {
  const std::size_t n = v.size();
  const std::size_t mid = n / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  const double upper = v[mid];
  if (n % 2 == 1) return upper;
  return 0.5 * (upper + *std::max_element(v.begin(), v.begin() + mid));
}

// Similarity weights of the training samples for one query. `w` (length n) is
// overwritten; samples outside the `k` nearest neighbours get a zero weight.
// With `robust_scale`, the scale of the distances is 1.4826 times their median
// absolute deviation instead of their standard deviation.
// Returns the number of non-zero weights.
arma::uword similarity_weights(const double* dist, const arma::uword n,
                               const arma::uword k, const double localization,
                               double* w, std::vector<arma::uword>& idx,
                               const bool robust_scale = false) {
  std::fill(w, w + n, 0.0);
  std::iota(idx.begin(), idx.end(), arma::uword(0));
  if (k < n) {
    std::nth_element(idx.begin(), idx.begin() + k, idx.end(),
                     [dist](arma::uword a, arma::uword b) {
                       return dist[a] < dist[b];
                     });
  }

  double d_min = std::numeric_limits<double>::infinity();
  double d_mean = 0.0;
  for (arma::uword s = 0; s < k; ++s) {
    const double d = dist[idx[s]];
    d_mean += d;
    d_min = std::min(d_min, d);
  }
  d_mean /= static_cast<double>(k);

  double ss = 0.0;
  for (arma::uword s = 0; s < k; ++s) {
    const double dev = dist[idx[s]] - d_mean;
    ss += dev * dev;
  }
  double d_sd = k > 1 ? std::sqrt(ss / static_cast<double>(k - 1)) : 0.0;

  if (robust_scale && k > 1) {
    std::vector<double> v(k);
    for (arma::uword s = 0; s < k; ++s) v[s] = dist[idx[s]];
    const double med = median_inplace(v);
    for (arma::uword s = 0; s < k; ++s) v[s] = std::abs(dist[idx[s]] - med);
    const double mad = 1.4826 * median_inplace(v);
    if (mad > 0.0) d_sd = mad;
  }

  // All distances equal: every sample is equally similar.
  if (!(d_sd > 0.0) || !std::isfinite(d_sd)) {
    for (arma::uword s = 0; s < k; ++s) w[idx[s]] = 1.0;
    return k;
  }

  // Shifting by the smallest distance rescales all weights by a constant,
  // which leaves the weighted PLS model unchanged but avoids underflow.
  const double bandwidth = d_sd * localization;
  arma::uword n_nonzero = 0;
  for (arma::uword s = 0; s < k; ++s) {
    const arma::uword i = idx[s];
    const double value = std::exp(-(dist[i] - d_min) / bandwidth);
    if (value >= weight_floor) {
      w[i] = value;
      ++n_nonzero;
    }
  }
  return n_nonzero;
}

}  // namespace

// Predictions of LW-PLS models with 1, ..., num_comp components.
//
// x, y     : training predictors (n x p) and outcomes (n x q), already
//            centred (and possibly scaled).
// new_x    : query predictors (m x p), transformed like x.
// dist_x   : coordinates of the training samples used for the distances
//            (n x r); x itself for Euclidean similarity.
// dist_new : the same coordinates for the queries (m x r).
// Returns an m x q x num_comp array of predictions on the scale of y.
// [[Rcpp::export(rng = false)]]
arma::cube lwpls_predict_cpp(const arma::mat& x, const arma::mat& y,
                             const arma::mat& new_x, const arma::mat& dist_x,
                             const arma::mat& dist_new, const int num_comp,
                             const double localization, const int neighbors,
                             const double tol) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword q = y.n_cols;
  const arma::uword m = new_x.n_rows;
  const arma::uword n_comp = static_cast<arma::uword>(std::max(num_comp, 1));
  const arma::uword k =
      std::min(n, static_cast<arma::uword>(std::max(neighbors, 1)));

  arma::cube out(m, q, n_comp, arma::fill::zeros);
  if (m == 0 || n == 0) return out;

  // Block size: keep the per-block working memory around 128 MB.
  const double per_query = 3.0 * n + 2.0 * p * n_comp + 1.0 * p * q + 4.0 * p +
                           1.0 * dist_x.n_cols;
  const arma::uword block = static_cast<arma::uword>(
      std::max(1.0, std::min(256.0, 16777216.0 / per_query)));

  // Products are written with the explicitly transposed training data: the
  // reference BLAS shipped with R is much slower for transposed products.
  const arma::mat xt = x.t();  // p x n
  const arma::mat yt = y.t();  // q x n
  const arma::vec x_sq = arma::sum(arma::square(x), 1);
  const arma::vec dist_sq = arma::sum(arma::square(dist_x), 1);
  std::vector<arma::uword> idx(n);
  arma::vec dist(n);

  for (arma::uword start = 0; start < m; start += block) {
    Rcpp::checkUserInterrupt();
    const arma::uword end = std::min(start + block, m);
    const arma::uword b = end - start;
    const arma::mat query = new_x.rows(start, end - 1).t();  // p x b

    // Distances and similarity weights -------------------------------------
    const arma::mat dist_query = dist_new.rows(start, end - 1).t();  // r x b
    arma::mat w = dist_x * dist_query;  // n x b; cross-products, overwritten
    const arma::rowvec query_sq = arma::sum(arma::square(dist_query), 0);
    std::vector<arma::uword> max_comp(b);
    for (arma::uword j = 0; j < b; ++j) {
      double* col = w.colptr(j);
      for (arma::uword i = 0; i < n; ++i) {
        const double d2 = dist_sq[i] + query_sq[j] - 2.0 * col[i];
        dist[i] = d2 > 0.0 ? std::sqrt(d2) : 0.0;
      }
      const arma::uword n_nonzero =
          similarity_weights(dist.memptr(), n, k, localization, col, idx);
      // A weighted-centred matrix with n_nonzero rows has rank < n_nonzero.
      max_comp[j] = std::min({n_comp, p, n_nonzero - 1});
    }
    const arma::rowvec w_sum = arma::sum(w, 0);

    // Local (weighted) means and centred queries ----------------------------
    arma::mat x_mean = xt * w;  // p x b
    x_mean.each_row() /= w_sum;
    arma::mat y_mean = yt * w;  // q x b
    y_mean.each_row() /= w_sum;
    const arma::mat query_c = query - x_mean;

    // Cross-products X_c' W Y_c, one p x q slice per query ------------------
    // Outcomes are centred explicitly, predictors through the identity
    // X_c' W Y_c = X' W Y_c - x_mean (1' W Y_c).
    arma::cube xy(p, q, b);
    arma::mat z(n, b);
    arma::rowvec ss_y(b, arma::fill::zeros);
    for (arma::uword c = 0; c < q; ++c) {
      for (arma::uword j = 0; j < b; ++j) {
        z.col(j) = w.col(j) % (y.col(c) - y_mean(c, j));
        ss_y[j] += arma::dot(z.col(j), y.col(c) - y_mean(c, j));
      }
      const arma::mat xz = xt * z;
      const arma::rowvec z_sum = arma::sum(z, 0);
      for (arma::uword j = 0; j < b; ++j) {
        xy.slice(j).col(c) = xz.col(j) - x_mean.col(j) * z_sum[j];
      }
    }

    // Total weighted sum of squares of the centred predictors (tolerances).
    arma::rowvec ss_x =
        (w.t() * x_sq).t() - w_sum % arma::sum(arma::square(x_mean), 0);
    ss_x.transform([](double v) { return v > 0.0 ? v : 0.0; });

    // Latent variables -------------------------------------------------------
    arma::mat pred = y_mean;  // q x b
    arma::cube r_hist(p, n_comp, b);
    arma::cube p_hist(p, n_comp, b);
    arma::vec xy_norm0(b);
    std::vector<char> active(b);
    for (arma::uword j = 0; j < b; ++j) {
      xy_norm0[j] = arma::norm(xy.slice(j), "fro");
      active[j] = max_comp[j] > 0 && std::isfinite(xy_norm0[j]) &&
                  xy_norm0[j] > cov_floor * std::sqrt(ss_x[j] * ss_y[j]);
    }

    arma::mat r_cur(p, b);
    arma::rowvec tt(b);
    for (arma::uword a = 0; a < n_comp; ++a) {
      // Weights w_a (dominant left singular vector of the deflated X'WY) and
      // the corresponding R weights r_a = w_a - R (P' w_a).
      r_cur.zeros();
      arma::uword n_active = 0;
      for (arma::uword j = 0; j < b; ++j) {
        if (!active[j]) continue;
        const arma::mat& xy_j = xy.slice(j);
        const double xy_norm = arma::norm(xy_j, "fro");
        if (a >= max_comp[j] || !(xy_norm > tol * xy_norm0[j])) {
          active[j] = 0;
          continue;
        }
        arma::vec w_a;
        if (q == 1) {
          w_a = xy_j.col(0) / xy_norm;
        } else {
          arma::vec eig_val;
          arma::mat eig_vec;
          if (!arma::eig_sym(eig_val, eig_vec, xy_j.t() * xy_j)) {
            // # nocov start
            // Defensive: the matrix is a small symmetric PSD matrix.
            active[j] = 0;
            continue;
          }  // # nocov end
          w_a = xy_j * eig_vec.col(q - 1);
          w_a /= arma::norm(w_a);
        }
        if (a > 0) {
          const arma::vec proj = p_hist.slice(j).head_cols(a).t() * w_a;
          r_cur.col(j) = w_a - r_hist.slice(j).head_cols(a) * proj;
        } else {
          r_cur.col(j) = w_a;
        }
        ++n_active;
      }

      if (n_active > 0) {
        // Scores of the training samples: t = X_c r = X r - (x_mean' r) 1.
        arma::mat t = x * r_cur;  // n x b
        for (arma::uword j = 0; j < b; ++j) {
          if (active[j]) {
            t.col(j) -= arma::dot(x_mean.col(j), r_cur.col(j));
            z.col(j) = w.col(j) % t.col(j);
            tt[j] = arma::dot(z.col(j), t.col(j));
            if (!(tt[j] > tol * ss_x[j]) || !std::isfinite(tt[j])) {
              // # nocov start
              // Defensive: exhausted predictors are normally caught by the
              // cross-product check above.
              active[j] = 0;
            }  // # nocov end
          }
          if (!active[j]) z.col(j).zeros();
        }

        // Loadings p_a = X_c' W t / (t' W t).
        const arma::mat xwt = xt * z;  // p x b
        for (arma::uword j = 0; j < b; ++j) {
          if (!active[j]) continue;
          const arma::vec& r_j = r_cur.col(j);
          const arma::vec x_load =
              (xwt.col(j) - x_mean.col(j) * arma::accu(z.col(j))) / tt[j];
          const arma::vec y_load = xy.slice(j).t() * r_j / tt[j];
          const double t_query = arma::dot(query_c.col(j), r_j);

          pred.col(j) += t_query * y_load;
          xy.slice(j) -= x_load * (tt[j] * y_load.t());
          r_hist.slice(j).col(a) = r_j;
          p_hist.slice(j).col(a) = x_load;
        }
      }

      // Queries whose local model is exhausted keep their last prediction.
      for (arma::uword j = 0; j < b; ++j) {
        for (arma::uword c = 0; c < q; ++c) {
          out(start + j, c, a) = pred(c, j);
        }
      }
    }
  }

  return out;
}

// -----------------------------------------------------------------------------
// Sparse and robust LW-PLS.
//
// Each query gets its own weighted NIPALS model, computed with explicit
// deflation of the training samples that have a non-zero similarity weight:
//
// * Sparsity (SNIPLS; Hoffmann et al. 2015): the weight vector of each
//   component is soft-thresholded at `sparsity * max|w|`, and the X loadings
//   are set to zero outside the variables selected so far.
// * Robustness (partial robust M-regression, PRM; Serneels et al. 2005): the
//   local model is refitted with case weights that are the products of the
//   similarity weights, residual weights and leverage weights (Fair or Hampel
//   functions), until the norm of the regression coefficients of the
//   normalized scores is stable. Centers and scales are similarity-weighted
//   medians and a similarity-weighted L1-median, so that the local model is
//   robust and local at the same time.

namespace {

enum RobustFun { ROBUST_NONE = 0, ROBUST_FAIR = 1, ROBUST_HAMPEL = 2 };

// Weighted median; equal weights give the usual median (the average of the
// two middle values when their number is even).
double weighted_median(const arma::vec& x, const arma::vec& w) {
  const arma::uword n = x.n_elem;
  const arma::uvec o = arma::sort_index(x);
  const double total = arma::accu(w);
  const double half = 0.5 * total;
  const double eps = 1e-9 * total;
  double cum = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    cum += w[o[i]];
    if (cum >= half - eps) {
      if (std::abs(cum - half) <= eps && i + 1 < n) {
        return 0.5 * (x[o[i]] + x[o[i + 1]]);
      }
      return x[o[i]];
    }
  }
  return x[o[n - 1]];  // # nocov (the cumulative weight reaches one half)
}

// Weighted quantile: the smallest value whose cumulative weight reaches the
// fraction `prob` of the total (the inverse of the weighted empirical
// distribution function; type 1 of stats::quantile() for equal weights).
double weighted_quantile(const arma::vec& x, const arma::vec& w, const double prob) {
  const arma::uword n = x.n_elem;
  const arma::uvec o = arma::sort_index(x);
  const double target = prob * arma::accu(w) * (1.0 - 1e-12);
  double cum = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    cum += w[o[i]];
    if (cum >= target) return x[o[i]];
  }
  return x[o[n - 1]];  // # nocov (the cumulative weight reaches the total)
}

// Weighted fraction of the values that are not larger than `value`.
double weighted_rank(const arma::vec& x, const arma::vec& w, const double value) {
  double below = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (x[i] <= value) below += w[i];
  }
  return below / arma::accu(w);
}

// Weighted median of non-negative values used as a scale: falls back to the
// weighted median of the positive values, then to one.
double weighted_scale(const arma::vec& x, const arma::vec& w) {
  const double s = weighted_median(x, w);
  if (s > 0.0) return s;
  const arma::uvec pos = arma::find(x > 0.0);
  if (pos.n_elem == 0) return 1.0;
  const arma::vec xp = x.elem(pos);
  const arma::vec wp = w.elem(pos);
  return weighted_median(xp, wp);
}

// Weighted L1-median (spatial median) of the rows of X, by Weiszfeld
// iterations started at the weighted mean. In one dimension it is the
// weighted median, where Weiszfeld iterations would converge slowly.
arma::vec weighted_l1median(const arma::mat& X, const arma::vec& w) {
  if (X.n_cols == 1) return arma::vec{weighted_median(X.col(0), w)};
  arma::vec m = X.t() * w / arma::accu(w);
  for (int it = 0; it < 1000; ++it) {
    const arma::mat diff = X.each_row() - m.t();
    arma::vec d = arma::sqrt(arma::sum(arma::square(diff), 1));
    const double d_max = d.max();
    if (!(d_max > 0.0)) break;
    d.transform([d_max](double v) { return std::max(v, 1e-12 * d_max); });
    const arma::vec v = w / d;
    const arma::vec m_new = X.t() * v / arma::accu(v);
    const double step = arma::norm(m_new - m);
    m = m_new;
    if (step <= 1e-12 * std::max(1.0, arma::norm(m))) break;
  }
  return m;
}

double fair_weight(const double z, const double c) {
  const double u = 1.0 + std::abs(z / c);
  return 1.0 / (u * u);
}

double hampel_weight(double z, const arma::vec& cut) {
  z = std::abs(z);
  if (z <= cut[0]) return 1.0;
  if (z <= cut[1]) return cut[0] / z;
  if (z <= cut[2]) return cut[0] * (cut[2] - z) / ((cut[2] - cut[1]) * z);
  return 0.0;
}

// Hampel cutoffs: quantiles of the chi distribution with `df` degrees of
// freedom at the probabilities `probs`.
arma::vec chi_cutoffs(const arma::vec& probs, const double df) {
  arma::vec cut(3);
  for (int i = 0; i < 3; ++i) cut[i] = std::sqrt(R::qchisq(probs[i], df, 1, 0));
  return cut;
}

// Weights of the rows of a residual matrix: distances of the robustly
// centered and scaled residuals, through the Fair or Hampel function.
arma::vec residual_weights(const arma::mat& E, const arma::vec& om,
                           const int fun, const double fair_c,
                           const arma::vec& probs, const double df) {
  const arma::uword k = E.n_rows;
  const arma::uword q = E.n_cols;
  arma::mat rc(E);
  for (arma::uword c = 0; c < q; ++c) {
    rc.col(c) -= weighted_median(E.col(c), om);
    double s = weighted_scale(arma::abs(rc.col(c)), om);
    if (fun == ROBUST_HAMPEL) s *= 1.4826;
    rc.col(c) /= s;
  }
  const arma::vec d = arma::sqrt(arma::sum(arma::square(rc), 1));
  arma::vec out(k);
  if (fun == ROBUST_FAIR) {
    const double md = weighted_scale(d, om);
    for (arma::uword i = 0; i < k; ++i) out[i] = fair_weight(d[i] / md, fair_c);
  } else {
    const arma::vec cut = chi_cutoffs(probs, std::max(df, 1.0));
    for (arma::uword i = 0; i < k; ++i) out[i] = hampel_weight(d[i], cut);
  }
  return out;
}

// Weights of the training samples from their distances in the score space.
arma::vec leverage_weights(const arma::mat& T, const arma::vec& tt,
                           const arma::vec& om, const int fun,
                           const double fair_c, const arma::vec& probs) {
  const arma::uword k = T.n_rows;
  const arma::uword a = T.n_cols;
  arma::vec out(k);
  if (fun == ROBUST_FAIR) {
    arma::mat Tn(T);
    for (arma::uword j = 0; j < a; ++j) Tn.col(j) /= std::sqrt(tt[j]);
    const arma::vec center = weighted_l1median(Tn, om);
    const arma::vec d =
        arma::sqrt(arma::sum(arma::square(Tn.each_row() - center.t()), 1));
    const double md = weighted_scale(d, om);
    for (arma::uword i = 0; i < k; ++i) out[i] = fair_weight(d[i] / md, fair_c);
  } else {
    arma::mat tc(T);
    for (arma::uword j = 0; j < a; ++j) {
      tc.col(j) -= weighted_median(T.col(j), om);
      tc.col(j) /= 1.4826 * weighted_scale(arma::abs(tc.col(j)), om);
    }
    const arma::vec d = arma::sqrt(arma::sum(arma::square(tc), 1));
    const arma::vec cut = chi_cutoffs(probs, static_cast<double>(a));
    for (arma::uword i = 0; i < k; ++i) out[i] = hampel_weight(d[i], cut);
  }
  return out;
}

struct NipalsFit {
  arma::mat T;       // k x ncomp training scores
  arma::mat C;       // q x ncomp Y loadings
  arma::vec tt;      // t' W t of each component
  arma::mat pred;    // q x a centered query predictions, 1..a components
  arma::uword ncomp = 0;  // number of components extracted
  // Diagnostics (only with `diagnostics = true`):
  arma::mat W;       // p x ncomp weight vectors
  arma::mat P;       // p x ncomp X loadings
  arma::vec t_query; // ncomp query scores
  arma::vec q_train; // k squared norms of the X residuals of the samples
  double q_query = 0.0;  // squared norm of the X residual of the query
};

// Weighted NIPALS on centered data, with optional SNIPLS sparsity. The query
// (centered like Xc) is projected along the way.
NipalsFit weighted_nipals(const arma::mat& Xc, const arma::mat& Yc,
                          const arma::vec& w, const arma::uword a,
                          const double sparsity, const arma::vec& xq,
                          const double tol, const bool diagnostics = false) {
  const arma::uword k = Xc.n_rows;
  const arma::uword p = Xc.n_cols;
  const arma::uword q = Yc.n_cols;
  NipalsFit fit;
  fit.T.zeros(k, a);
  fit.C.zeros(q, a);
  fit.tt.zeros(a);
  fit.pred.zeros(q, a);
  fit.ncomp = 0;
  if (diagnostics) {
    fit.W.zeros(p, a);
    fit.P.zeros(p, a);
    fit.t_query.zeros(a);
  }

  arma::mat Xa(Xc);
  arma::mat Ya(Yc);
  arma::vec xqa(xq);
  arma::vec cur(q, arma::fill::zeros);
  std::vector<char> active(p, 0);

  const double ss_x = arma::accu(arma::sum(arma::square(Xc), 1) % w);
  const double ss_y = arma::accu(arma::sum(arma::square(Yc), 1) % w);
  double m_norm0 = 0.0;

  for (arma::uword j = 0; j < a; ++j) {
    const arma::mat M = Xa.t() * (Ya.each_col() % w);  // p x q
    const double m_norm = arma::norm(M, "fro");
    if (j == 0) {
      m_norm0 = m_norm;
      if (!std::isfinite(m_norm) ||
          !(m_norm > cov_floor * std::sqrt(ss_x * ss_y))) {
        break;
      }
    } else if (!(m_norm > tol * m_norm0)) {
      break;
    }

    arma::vec wv;
    if (q == 1) {
      wv = M.col(0) / m_norm;
    } else {
      arma::vec eig_val;
      arma::mat eig_vec;
      if (!arma::eig_sym(eig_val, eig_vec, M.t() * M)) break;  // # nocov
      wv = M * eig_vec.col(q - 1);
      wv /= arma::norm(wv);
    }

    if (sparsity > 0.0) {
      const double thr = sparsity * arma::abs(wv).max();
      for (arma::uword i = 0; i < p; ++i) {
        const double aw = std::abs(wv[i]);
        if (aw >= thr) active[i] = 1;
        wv[i] = aw > thr ? (wv[i] > 0.0 ? aw - thr : thr - aw) : 0.0;
      }
      const double nrm = arma::norm(wv);
      if (!(nrm > 0.0)) break;
      wv /= nrm;
    }

    const arma::vec t = Xa * wv;
    const double tt = arma::dot(w % t, t);
    if (!(tt > tol * ss_x) || !std::isfinite(tt)) break;
    arma::vec pv = Xa.t() * (w % t) / tt;
    if (sparsity > 0.0) {
      for (arma::uword i = 0; i < p; ++i) {
        if (!active[i]) pv[i] = 0.0;
      }
    }
    const arma::vec cv = Ya.t() * (w % t) / tt;

    const double tq = arma::dot(xqa, wv);
    cur += tq * cv;
    xqa -= tq * pv;
    Xa -= t * pv.t();
    Ya -= t * cv.t();

    fit.T.col(j) = t;
    fit.C.col(j) = cv;
    fit.tt[j] = tt;
    fit.pred.col(j) = cur;
    if (diagnostics) {
      fit.W.col(j) = wv;
      fit.P.col(j) = pv;
      fit.t_query[j] = tq;
    }
    fit.ncomp = j + 1;
  }
  for (arma::uword j = fit.ncomp; j < a; ++j) fit.pred.col(j) = cur;
  fit.T = fit.T.head_cols(fit.ncomp);
  fit.C = fit.C.head_cols(fit.ncomp);
  fit.tt = fit.tt.head(fit.ncomp);
  if (diagnostics) {
    fit.W = fit.W.head_cols(fit.ncomp);
    fit.P = fit.P.head_cols(fit.ncomp);
    fit.t_query = fit.t_query.head(fit.ncomp);
    fit.q_train = arma::sum(arma::square(Xa), 1);
    fit.q_query = arma::dot(xqa, xqa);
  }
  return fit;
}

// The training samples, centers and initial weights of one query.
struct LocalProblem {
  arma::uvec rows;     // training samples with a non-zero similarity weight
  arma::vec om;        // their similarity weights
  arma::mat Xc;        // their centered predictors
  arma::mat Yc;        // their centered outcomes
  arma::rowvec cy;     // center of the outcomes
  arma::vec xqc;       // centered query
  arma::vec w_init;    // initial robust weights (robust models)
  arma::uword cap = 0; // maximum number of components
  double nearest = 0.0;
};

struct Settings {
  arma::uword n_comp;
  arma::uword k_nn;
  double localization;
  double sparsity;
  int robust;
  double fair_c;
  arma::vec hampel_probs;
  int max_iter;
  bool classification;
  double df_res;
  double tol;
};

LocalProblem prepare_local(const arma::mat& x, const arma::mat& y,
                           const arma::vec& dist, const arma::vec& xq,
                           const Settings& set, std::vector<arma::uword>& idx,
                           arma::vec& omega) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword q = y.n_cols;
  LocalProblem lp;
  similarity_weights(dist.memptr(), n, set.k_nn, set.localization,
                     omega.memptr(), idx, set.robust != ROBUST_NONE);
  lp.nearest = dist.min();
  lp.rows = arma::find(omega > 0.0);
  lp.om = omega.elem(lp.rows);
  const arma::mat Xs = x.rows(lp.rows);
  const arma::mat Ys = y.rows(lp.rows);
  const arma::uword k = lp.rows.n_elem;
  // A weighted-centred matrix with k rows has rank < k.
  lp.cap = std::min({set.n_comp, p, k > 0 ? k - 1 : 0});

  if (set.robust == ROBUST_NONE) {
    const arma::rowvec cx = (Xs.t() * lp.om).t() / arma::accu(lp.om);
    lp.cy = (Ys.t() * lp.om).t() / arma::accu(lp.om);
    lp.Xc = Xs.each_row() - cx;
    lp.Yc = Ys.each_row() - lp.cy;
    lp.xqc = xq - cx.t();
    return lp;
  }

  // Robust local model (PRM): L1-median and medians weighted by similarity.
  const arma::vec cx = weighted_l1median(Xs, lp.om);
  lp.cy.set_size(q);
  for (arma::uword c = 0; c < q; ++c) {
    lp.cy[c] = set.classification ? arma::dot(Ys.col(c), lp.om) / arma::accu(lp.om)
                                  : weighted_median(Ys.col(c), lp.om);
  }
  lp.Xc = Xs.each_row() - cx.t();
  lp.Yc = Ys.each_row() - lp.cy;
  lp.xqc = xq - cx;

  // Initial weights from the distances to the centers.
  const arma::vec dx = arma::sqrt(arma::sum(arma::square(lp.Xc), 1));
  const double mdx = weighted_scale(dx, lp.om);
  lp.w_init.set_size(k);
  if (set.robust == ROBUST_FAIR) {
    for (arma::uword i = 0; i < k; ++i) lp.w_init[i] = fair_weight(dx[i] / mdx, set.fair_c);
  } else {
    arma::vec cut(3);
    for (int i = 0; i < 3; ++i) cut[i] = R::qnorm(set.hampel_probs[i], 0.0, 1.0, 1, 0);
    for (arma::uword i = 0; i < k; ++i) lp.w_init[i] = hampel_weight(dx[i] / mdx, cut);
  }
  if (!set.classification) {
    lp.w_init %= residual_weights(lp.Yc, lp.om, set.robust, set.fair_c,
                                  set.hampel_probs, set.df_res);
  }
  lp.w_init.transform([](double v) { return std::max(v, 1e-6); });
  return lp;
}

struct LocalFit {
  NipalsFit fit;
  arma::vec w_fit;     // case weights of the final fit
  arma::rowvec pred;   // prediction of the query (scale of y)
};

// Local model of one query with `a` components (1 <= a <= cap).
LocalFit fit_local(const LocalProblem& lp, const arma::uword a,
                   const Settings& set, const bool diagnostics) {
  const double prm_tol = 0.01;  // relative change of the coefficient norm
  const arma::uword q = lp.Yc.n_cols;
  LocalFit out;
  if (set.robust == ROBUST_NONE) {
    out.fit = weighted_nipals(lp.Xc, lp.Yc, lp.om, a, set.sparsity, lp.xqc,
                              set.tol, diagnostics);
    out.w_fit = lp.om;
    out.pred = lp.cy + out.fit.pred.col(a - 1).t();
    return out;
  }

  arma::vec wr(lp.w_init);
  double gamma = 1e5;
  double diff = 1.0;
  int iter = 1;
  while (diff > prm_tol && iter <= set.max_iter) {
    out.w_fit = lp.om % wr;
    out.fit = weighted_nipals(lp.Xc, lp.Yc, out.w_fit, a, set.sparsity,
                              lp.xqc, set.tol, diagnostics);
    double g = 0.0;
    for (arma::uword h = 0; h < out.fit.ncomp; ++h) {
      g += out.fit.tt[h] * arma::dot(out.fit.C.col(h), out.fit.C.col(h));
    }
    g = std::sqrt(g);
    diff = g > 0.0 ? std::abs(g - gamma) / g : 0.0;
    gamma = g;
    if (out.fit.ncomp == 0) break;
    const arma::mat E = lp.Yc - out.fit.T * out.fit.C.t();
    wr = residual_weights(E, lp.om, set.robust, set.fair_c, set.hampel_probs, set.df_res) %
         leverage_weights(out.fit.T, out.fit.tt, lp.om, set.robust, set.fair_c,
                          set.hampel_probs);
    wr.transform([](double v) { return std::max(v, 1e-6); });
    ++iter;
  }
  // Intercept from the residuals of the last fit: weighted medians, or
  // weighted means for class indicators so that predictions sum to one.
  const arma::mat E = out.fit.ncomp > 0
                          ? arma::mat(lp.Yc - out.fit.T * out.fit.C.t())
                          : lp.Yc;
  out.pred.set_size(q);
  for (arma::uword c = 0; c < q; ++c) {
    const double adj = set.classification
                           ? arma::dot(E.col(c), out.w_fit) / arma::accu(out.w_fit)
                           : weighted_median(E.col(c), lp.om);
    out.pred[c] = lp.cy[c] + adj + out.fit.pred(c, a - 1);
  }
  return out;
}

Settings make_settings(const arma::uword n, const arma::uword q,
                       const arma::uword n_comp, const double localization,
                       const int neighbors, const double sparsity,
                       const int robust, const double fair_c,
                       const arma::vec& hampel_probs, const int max_iter,
                       const bool classification, const double tol) {
  Settings set;
  set.n_comp = n_comp;
  set.k_nn = std::min(n, static_cast<arma::uword>(std::max(neighbors, 1)));
  set.localization = localization;
  set.sparsity = sparsity;
  set.robust = robust;
  set.fair_c = fair_c;
  set.hampel_probs = hampel_probs;
  set.max_iter = max_iter;
  set.classification = classification;
  set.df_res = classification ? static_cast<double>(q) - 1.0 : static_cast<double>(q);
  set.tol = tol;
  return set;
}

// Distances between the query `j` and the training samples.
void query_distances(const arma::mat& dist_x, const arma::vec& dist_sq,
                     const arma::mat& dist_new, const arma::uword j,
                     arma::vec& dist) {
  const arma::vec dq = dist_new.row(j).t();
  const arma::vec cross = dist_x * dq;
  const double dq_sq = arma::dot(dq, dq);
  for (arma::uword i = 0; i < dist_x.n_rows; ++i) {
    const double d2 = dist_sq[i] + dq_sq - 2.0 * cross[i];
    dist[i] = d2 > 0.0 ? std::sqrt(d2) : 0.0;
  }
}

}  // namespace

// Predictions of sparse and/or robust LW-PLS models.
//
// Arguments as for lwpls_predict_cpp(), and:
// comps          : numbers of components to predict (1-based). Non-robust
//                  models give all numbers up to max(comps) from one fit;
//                  robust models are refitted for each element of `comps`.
// sparsity       : SNIPLS threshold (0 for none).
// robust         : 0 (none), 1 (Fair) or 2 (Hampel).
// fair_c         : constant of the Fair function.
// hampel_probs   : probabilities of the three Hampel cutoffs.
// max_iter       : maximum number of PRM fits (at least one).
// classification : whether y holds class indicators.
// Returns an m x q x max(comps) array; slices not requested are NaN.
// [[Rcpp::export(rng = false)]]
arma::cube lwpls_general_cpp(const arma::mat& x, const arma::mat& y,
                             const arma::mat& new_x, const arma::mat& dist_x,
                             const arma::mat& dist_new,
                             const arma::uvec& comps, const double localization,
                             const int neighbors, const double sparsity,
                             const int robust, const double fair_c,
                             const arma::vec& hampel_probs, const int max_iter,
                             const bool classification, const double tol) {
  const arma::uword n = x.n_rows;
  const arma::uword q = y.n_cols;
  const arma::uword m = new_x.n_rows;
  const arma::uword a_max = comps.max();

  arma::cube out(m, q, a_max);
  out.fill(arma::datum::nan);
  if (m == 0 || n == 0) return out;

  const Settings set = make_settings(n, q, a_max, localization, neighbors, sparsity,
                                     robust, fair_c, hampel_probs, max_iter,
                                     classification, tol);
  const arma::vec dist_sq = arma::sum(arma::square(dist_x), 1);
  std::vector<arma::uword> idx(n);
  arma::vec dist(n);
  arma::vec omega(n);

  for (arma::uword j = 0; j < m; ++j) {
    if (j % 16 == 0) Rcpp::checkUserInterrupt();
    query_distances(dist_x, dist_sq, dist_new, j, dist);
    const LocalProblem lp = prepare_local(x, y, dist, new_x.row(j).t(), set, idx, omega);

    if (robust == ROBUST_NONE) {
      // One fit gives all the numbers of components.
      arma::mat pred(q, a_max);
      if (lp.cap > 0) {
        const LocalFit lf = fit_local(lp, lp.cap, set, false);
        for (arma::uword a = 0; a < a_max; ++a) {
          pred.col(a) = lp.cy.t() + lf.fit.pred.col(std::min(a, lp.cap - 1));
        }
      } else {
        pred.each_col() = lp.cy.t();
      }
      for (arma::uword a = 0; a < a_max; ++a) {
        for (arma::uword c = 0; c < q; ++c) out(j, c, a) = pred(c, a);
      }
      continue;
    }

    for (arma::uword ci = 0; ci < comps.n_elem; ++ci) {
      const arma::uword a_req = comps[ci];
      const arma::uword a = std::min(a_req, lp.cap);
      const arma::rowvec pred = a > 0 ? fit_local(lp, a, set, false).pred : lp.cy;
      for (arma::uword c = 0; c < q; ++c) out(j, c, a_req - 1) = pred[c];
    }
  }
  return out;
}

// Local models of the queries, for diagnostics.
//
// Arguments as for lwpls_general_cpp(), with a single number of components,
// and the coverage `level` of the reliability limits.
// Returns a list with, for each query, the prediction, the regression
// coefficients (on the standardized scale), the VIP, the selected predictors,
// the reliability statistics (columns: effective number of neighbors,
// distance to the nearest training sample, number of components, T2 of the
// query and its theoretical limit, Q of the query and its theoretical limit,
// the empirical limits of T2 and Q, and the percentile ranks of the T2 and Q
// of the query among the training samples of its local model), and, over
// all queries, the sums of the similarity
// weights and of the similarity weights times the robust weights of the
// training samples, the number of queries that select each predictor in each
// component, and the number of queries whose local model has each component.
// [[Rcpp::export(rng = false)]]
Rcpp::List lwpls_local_cpp(const arma::mat& x, const arma::mat& y,
                           const arma::mat& new_x, const arma::mat& dist_x,
                           const arma::mat& dist_new, const int num_comp,
                           const double localization, const int neighbors,
                           const double sparsity, const int robust,
                           const double fair_c, const arma::vec& hampel_probs,
                           const int max_iter, const bool classification,
                           const double tol, const double level) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword q = y.n_cols;
  const arma::uword m = new_x.n_rows;
  const Settings set = make_settings(n, q, static_cast<arma::uword>(std::max(num_comp, 1)),
                                     localization, neighbors, sparsity, robust,
                                     fair_c, hampel_probs, max_iter,
                                     classification, tol);

  arma::mat pred(m, q, arma::fill::zeros);
  arma::cube coef(m, p, q, arma::fill::zeros);
  arma::mat vip(m, p);
  vip.fill(arma::datum::nan);
  arma::mat selected(m, p, arma::fill::zeros);
  arma::mat reliability(m, 11);
  reliability.fill(arma::datum::nan);
  arma::vec omega_sum(n, arma::fill::zeros);
  arma::vec omega_robust_sum(n, arma::fill::zeros);
  arma::mat selected_comp(p, set.n_comp, arma::fill::zeros);
  arma::vec comp_count(set.n_comp, arma::fill::zeros);

  const arma::vec dist_sq = arma::sum(arma::square(dist_x), 1);
  std::vector<arma::uword> idx(n);
  arma::vec dist(n);
  arma::vec omega(n);

  for (arma::uword j = 0; j < m; ++j) {
    if (j % 16 == 0) Rcpp::checkUserInterrupt();
    query_distances(dist_x, dist_sq, dist_new, j, dist);
    const LocalProblem lp = prepare_local(x, y, dist, new_x.row(j).t(), set, idx, omega);
    const double om_sum = arma::accu(lp.om);
    reliability(j, 0) = om_sum * om_sum / arma::dot(lp.om, lp.om);
    reliability(j, 1) = lp.nearest;

    const arma::uword a = std::min(set.n_comp, lp.cap);
    reliability(j, 2) = static_cast<double>(a);
    if (a == 0) {
      pred.row(j) = lp.cy;
      omega_sum.elem(lp.rows) += lp.om;
      omega_robust_sum.elem(lp.rows) += lp.om;
      continue;
    }
    const LocalFit lf = fit_local(lp, a, set, true);
    const NipalsFit& fit = lf.fit;
    pred.row(j) = lf.pred;
    omega_sum.elem(lp.rows) += lp.om;
    omega_robust_sum.elem(lp.rows) += lf.w_fit;
    if (fit.ncomp == 0) continue;
    const arma::uword h = fit.ncomp;

    // Regression coefficients B = W (P'W)^-1 C' and VIP.
    arma::mat ptw = fit.P.t() * fit.W;
    arma::mat inv;
    if (!arma::inv(inv, ptw)) inv = arma::pinv(ptw);  // # nocov
    const arma::mat B = fit.W * inv * fit.C.t();       // p x q
    for (arma::uword c = 0; c < q; ++c) {
      for (arma::uword i = 0; i < p; ++i) coef(j, i, c) = B(i, c);
    }
    arma::vec ssy(h);
    for (arma::uword a2 = 0; a2 < h; ++a2) {
      ssy[a2] = fit.tt[a2] * arma::dot(fit.C.col(a2), fit.C.col(a2));
    }
    const double ssy_tot = arma::accu(ssy);
    for (arma::uword a2 = 0; a2 < h; ++a2) comp_count[a2] += 1.0;
    for (arma::uword i = 0; i < p; ++i) {
      double num = 0.0;
      bool sel = false;
      for (arma::uword a2 = 0; a2 < h; ++a2) {
        num += ssy[a2] * fit.W(i, a2) * fit.W(i, a2);
        if (fit.W(i, a2) != 0.0) {
          sel = true;
          selected_comp(i, a2) += 1.0;
        }
      }
      vip(j, i) = ssy_tot > 0.0 ? std::sqrt(static_cast<double>(p) * num / ssy_tot) : 0.0;
      selected(j, i) = sel ? 1.0 : 0.0;
    }

    // Hotelling T2 and Q of the query and of the training samples in the local
    // model, with the weights of the fit.
    const arma::vec& w = lf.w_fit;
    const double w_sum = arma::accu(w);
    const double n_eff = w_sum * w_sum / arma::dot(w, w);
    double t2 = 0.0;
    arma::vec t2_train(fit.T.n_rows, arma::fill::zeros);
    for (arma::uword a2 = 0; a2 < h; ++a2) {
      const arma::vec t = fit.T.col(a2);
      const double mean = arma::dot(w, t) / w_sum;
      double var = arma::dot(w, arma::square(t - mean)) / w_sum;
      if (n_eff > 1.0) var *= n_eff / (n_eff - 1.0);
      if (!(var > 0.0)) continue;
      const double dev = fit.t_query[a2] - mean;
      t2 += dev * dev / var;
      t2_train += arma::square(t - mean) / var;
    }
    reliability(j, 3) = t2;
    reliability(j, 5) = fit.q_query;

    // Theoretical limits, from the effective number of samples of the fit.
    const double hd = static_cast<double>(h);
    if (n_eff > hd + 1.0) {
      reliability(j, 4) = hd * (n_eff * n_eff - 1.0) / (n_eff * (n_eff - hd)) *
                          R::qf(level, hd, n_eff - hd, 1, 0);
    }
    const double q_mean = arma::dot(w, fit.q_train) / w_sum;
    const double q_var = arma::dot(w, arma::square(fit.q_train - q_mean)) / w_sum;
    if (q_mean > 0.0 && q_var > 0.0) {
      // Box (1954) approximation, as in Nomikos & MacGregor (1995).
      const double g = q_var / (2.0 * q_mean);
      const double df = 2.0 * q_mean * q_mean / q_var;
      reliability(j, 6) = g * R::qchisq(level, df, 1, 0);
    } else if (q_mean > 0.0) {
      // # nocov start
      // Defensive: all the residuals are equal, a degenerate distribution.
      reliability(j, 6) = q_mean;
    }  // # nocov end

    // Empirical limits and percentile ranks among the training samples of
    // the local model.
    reliability(j, 7) = weighted_quantile(t2_train, w, level);
    reliability(j, 8) = weighted_quantile(fit.q_train, w, level);
    reliability(j, 9) = weighted_rank(t2_train, w, t2);
    reliability(j, 10) = weighted_rank(fit.q_train, w, fit.q_query);
  }

  return Rcpp::List::create(
      Rcpp::Named("pred") = pred,
      Rcpp::Named("coef") = coef,
      Rcpp::Named("vip") = vip,
      Rcpp::Named("selected") = selected,
      Rcpp::Named("reliability") = reliability,
      Rcpp::Named("omega_sum") = omega_sum,
      Rcpp::Named("omega_robust_sum") = omega_robust_sum,
      Rcpp::Named("selected_comp") = selected_comp,
      Rcpp::Named("comp_count") = comp_count);
}
