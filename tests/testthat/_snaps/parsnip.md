# model specification

    Code
      spec
    Output
      Locally-Weighted PLS Model Specification (unknown mode)
      
      Main Arguments:
        num_comp = 3
        localization = 0.5
      
      Computational engine: lwpls 
      

---

    Code
      parsnip::translate(parsnip::set_engine(parsnip::set_mode(lwpls(num_comp = 3,
        neighbors = 50), "classification"), "lwpls", scale = FALSE))
    Output
      Locally-Weighted PLS Model Specification (classification)
      
      Main Arguments:
        num_comp = 3
        neighbors = 50
      
      Engine-Specific Arguments:
        scale = FALSE
      
      Computational engine: lwpls 
      
      Model fit template:
      lwpls::lwpls_fit(x = missing_arg(), y = missing_arg(), num_comp = 3, 
          neighbors = 50, scale = FALSE)

# update method

    Code
      update(spec, localization = 0.25)
    Output
      Locally-Weighted PLS Model Specification (unknown mode)
      
      Main Arguments:
        num_comp = 3
        localization = 0.25
      
      Computational engine: lwpls 
      

---

    Code
      update(spec, localization = 0.25, fresh = TRUE)
    Output
      Locally-Weighted PLS Model Specification (unknown mode)
      
      Main Arguments:
        localization = 0.25
      
      Computational engine: lwpls 
      

---

    Code
      update(spec, param, scale = FALSE)
    Condition
      Error in `update_dot_check()`:
      ! The extra argument `scale` will be ignored.

# argument checks happen at fit time

    Code
      parsnip::fit(parsnip::set_mode(lwpls(localization = -1), "regression"), y ~ .,
      data = train)
    Condition
      Error in `parsnip::fit()`:
      ! `localization` must be a single positive number, not -1.

---

    Code
      parsnip::fit(parsnip::set_mode(lwpls(num_comp = 0), "regression"), y ~ ., data = train)
    Condition
      Error in `parsnip::fit()`:
      ! `num_comp` must be a single whole number >= 1, not 0.

# the similarity index is passed on

    Code
      parsnip::translate(spec)
    Output
      Locally-Weighted PLS Model Specification (regression)
      
      Main Arguments:
        num_comp = 2
        similarity = covariance
      
      Computational engine: lwpls 
      
      Model fit template:
      lwpls::lwpls_fit(x = missing_arg(), y = missing_arg(), num_comp = 2, 
          similarity = "covariance")

---

    Code
      parsnip::fit(parsnip::set_mode(lwpls(similarity = "cosine"), "regression"), y ~
        ., data = train)
    Condition
      Error in `parsnip::fit()`:
      ! `similarity` must be one of "euclidean" or "covariance", not "cosine".

# multi_predict() for regression

    Code
      parsnip::multi_predict(fit, newdata = test)
    Condition
      Error in `parsnip::multi_predict()`:
      ! Please use `new_data` instead of `newdata`.

---

    Code
      parsnip::multi_predict(fit, test, type = "prob")
    Condition
      Error in `parsnip::multi_predict()`:
      ! `type` must be one of "numeric", not "prob".

---

    Code
      parsnip::multi_predict(fit, test, num_comp = -1)
    Condition
      Error in `parsnip::multi_predict()`:
      ! `num_comp` must be a vector of whole numbers >= 1.

# classification fit and predictions

    Code
      predict(fit, new, type = "numeric")
    Condition
      Error in `predict()`:
      ! For numeric predictions, the object should be a regression model.

