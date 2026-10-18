# print method for sparse and robust models

    Code
      lwpls_fit(train[x_names], train$y, sparsity = 0.3, robust = TRUE)
    Output
      Locally-weighted PLS (regression)
      
      Training samples: 65 
      Predictors:       6 (standardized) 
      Components:       2 
      Localization:     1 
      Similarity:       Euclidean 
      Neighbors:        all (65) 
      Sparsity:         0.3 (SNIPLS)
      Robust:           PRM, Fair weights (c = 4)

---

    Code
      lwpls_fit(train[x_names], train$y, robust = TRUE, weight_function = "hampel")
    Output
      Locally-weighted PLS (regression)
      
      Training samples: 65 
      Predictors:       6 (standardized) 
      Components:       2 
      Localization:     1 
      Similarity:       Euclidean 
      Neighbors:        all (65) 
      Robust:           PRM, Hampel weights (p = 0.950, 0.975, 0.999)

# bad sparse and robust settings are rejected

    Code
      lwpls_fit(x, y, sparsity = 1)
    Condition
      Error in `lwpls_fit()`:
      ! `sparsity` must be a single number in [0, 1), not 1.

---

    Code
      lwpls_fit(x, y, sparsity = -0.1)
    Condition
      Error in `lwpls_fit()`:
      ! `sparsity` must be a single number in [0, 1), not -0.1.

---

    Code
      lwpls_fit(x, y, robust = "yes")
    Condition
      Error in `lwpls_fit()`:
      ! `robust` must be `TRUE` or `FALSE`, not a string.

---

    Code
      lwpls_fit(x, y, robust = TRUE, weight_function = "huber")
    Condition
      Error in `lwpls_fit()`:
      ! `weight_function` must be one of "fair" or "hampel", not "huber".

---

    Code
      lwpls_fit(x, y, robust = TRUE, robust_constant = -1)
    Condition
      Error in `lwpls_fit()`:
      ! `robust_constant` must be a single positive number, not -1.

---

    Code
      lwpls_fit(x, y, robust = TRUE, weight_function = "hampel", robust_constant = c(
        0.9, 0.8, 0.99))
    Condition
      Error in `lwpls_fit()`:
      ! `robust_constant` must be three increasing probabilities for the Hampel function, not a double vector.
      i For example, `c(0.95, 0.975, 0.999)`.

---

    Code
      lwpls_fit(x, y, robust = TRUE, max_iter = 0)
    Condition
      Error in `lwpls_fit()`:
      ! `max_iter` must be a single whole number >= 1, not 0.

# parsnip passes the sparse and robust settings

    Code
      parsnip::translate(spec)
    Output
      Locally-Weighted PLS Model Specification (regression)
      
      Main Arguments:
        num_comp = 3
        sparsity = 0.3
      
      Engine-Specific Arguments:
        robust = TRUE
        weight_function = hampel
      
      Computational engine: lwpls 
      
      Model fit template:
      lwpls::lwpls_fit(x = missing_arg(), y = missing_arg(), num_comp = 3, 
          sparsity = 0.3, robust = TRUE, weight_function = "hampel")

---

    Code
      parsnip::fit(parsnip::set_mode(lwpls(sparsity = 2), "regression"), y ~ ., data = train[
        c(x_names, "y")])
    Condition
      Error in `parsnip::fit()`:
      ! `sparsity` must be a single number in [0, 1), not 2.

# sparsity parameter

    Code
      sparsity()
    Message
      Sparsity Threshold (quantitative)
      Range: [0, 0.9]

# trimmed RMSE

    Code
      trimmed_rmse_vec(1:3, 1:3, trim = 1)
    Condition
      Error in `trimmed_rmse_vec()`:
      ! `trim` must be a single number in [0, 1), not 1.

