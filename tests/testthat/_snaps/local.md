# missing values and other interfaces

    Code
      lwpls_local(list())
    Condition
      Error in `lwpls_local()`:
      ! `lwpls_local()` is not defined for an empty list.

---

    Code
      lwpls_local(fit, test[x_names], foo = 1)
    Condition
      Error in `lwpls_local()`:
      ! `...` must be empty.
      x Problematic argument:
      * foo = 1

---

    Code
      lwpls_local(fit, test[x_names], level = 1)
    Condition
      Error in `lwpls_local()`:
      ! `level` must be a single number in (0, 1), not 1.

---

    Code
      lwpls_local(fit, test[x_names], limits = "exact")
    Condition
      Error in `lwpls_local()`:
      ! `limits` must be one of "empirical" or "theoretical", not "exact".

---

    Code
      lwpls_local(lm_fit, test)
    Condition
      Error in `lwpls_local()`:
      ! `object` must be a fitted `lwpls()` model.

# spectral axis

    Code
      spectral_axis(c("a", "b"), wavelength = 1:3)
    Condition
      Error:
      ! `wavelength` must be 2 unique finite numbers, one per predictor.

# print method

    Code
      lwpls_local(fit, test[x_names])
    Output
      LW-PLS local models
      
      Queries:           15 
      Components:        2 
      Reliability:       12 inside, 3 high Q (empirical limits, 95%) 
      Selected:          5 of 6 predictors (median)
      Robust weights:    median 0.33 | lowest: #62 (0.032), #56 (0.041), #28 (0.046) 

# plots

    Code
      plot_reliability(loc, reference = "percentile", log = TRUE)
    Condition
      Error in `plot_reliability()`:
      ! `log = TRUE` is not available for percentile ranks.

# plots for several outcomes and classes

    Code
      plot_coefficients(loc, outcome = "z")
    Condition
      Error in `plot_coefficients()`:
      ! `outcome` must be one of "y" or "y2", not "z".

# plots reject unsuitable models

    Code
      plot_robust_weights(loc)
    Condition
      Error in `plot_robust_weights()`:
      ! Robust weights are only available for robust models.
      i Fit the model with `robust = TRUE`.

---

    Code
      plot_selection(loc)
    Condition
      Error in `plot_selection()`:
      ! Predictor selection is only available for sparse models.
      i Fit the model with `sparsity > 0`.

---

    Code
      plot_reliability(list())
    Condition
      Error in `plot_reliability()`:
      ! `object` must be the result of `lwpls_local()`, not an empty list.

