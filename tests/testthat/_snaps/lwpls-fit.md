# wide data (p > n) is supported

    Code
      fit <- lwpls_fit(wide[1:15, ], y[1:15], num_comp = 20)
    Condition
      Warning:
      ! `num_comp` = 20 is larger than the maximum number of components (14).
      i 14 components will be used.

# print method

    Code
      lwpls_fit(train[x_names], train$y, num_comp = 2, neighbors = 50)
    Output
      Locally-weighted PLS (regression)
      
      Training samples: 90 
      Predictors:       6 (standardized) 
      Components:       2 
      Localization:     1 
      Similarity:       Euclidean 
      Neighbors:        50 

---

    Code
      lwpls_fit(Species ~ ., data = iris, scale = FALSE)
    Output
      Locally-weighted PLS (classification)
      
      Training samples: 150 
      Predictors:       4 
      Classes:          setosa, versicolor, virginica 
      Components:       2 
      Localization:     1 
      Similarity:       Euclidean 
      Neighbors:        all (150) 

---

    Code
      lwpls_fit(y + y2 ~ ., data = train, num_comp = 2)
    Output
      Locally-weighted PLS (regression)
      
      Training samples: 90 
      Predictors:       6 (standardized) 
      Outcomes:         y, y2 
      Components:       2 
      Localization:     1 
      Similarity:       Euclidean 
      Neighbors:        all (90) 

---

    Code
      lwpls_fit(y ~ ., data = train, similarity = "covariance")
    Output
      Locally-weighted PLS (regression)
      
      Training samples: 90 
      Predictors:       7 (standardized) 
      Outcomes:         y 
      Components:       2 
      Localization:     1 
      Similarity:       covariance-based (CbLW-PLS) 
      Neighbors:        all (90) 

# bad hyperparameters are rejected

    Code
      lwpls_fit(x, y, num_comp = 0)
    Condition
      Error in `lwpls_fit()`:
      ! `num_comp` must be a single whole number >= 1, not 0.

---

    Code
      lwpls_fit(x, y, num_comp = 1.5)
    Condition
      Error in `lwpls_fit()`:
      ! `num_comp` must be a single whole number >= 1, not 1.5.

---

    Code
      lwpls_fit(x, y, num_comp = 1:2)
    Condition
      Error in `lwpls_fit()`:
      ! `num_comp` must be a single whole number >= 1, not an integer vector.

---

    Code
      lwpls_fit(x, y, localization = -1)
    Condition
      Error in `lwpls_fit()`:
      ! `localization` must be a single positive number, not -1.

---

    Code
      lwpls_fit(x, y, localization = Inf)
    Condition
      Error in `lwpls_fit()`:
      ! `localization` must be a single positive number, not Inf.

---

    Code
      lwpls_fit(x, y, localization = "small")
    Condition
      Error in `lwpls_fit()`:
      ! `localization` must be a single positive number, not "small".

---

    Code
      lwpls_fit(x, y, neighbors = 0)
    Condition
      Error in `lwpls_fit()`:
      ! `neighbors` must be a single whole number >= 1, not 0.

---

    Code
      lwpls_fit(x, y, similarity = "mahalanobis")
    Condition
      Error in `lwpls_fit()`:
      ! `similarity` must be one of "euclidean" or "covariance", not "mahalanobis".

---

    Code
      lwpls_fit(x, y, scale = "yes")
    Condition
      Error in `lwpls_fit()`:
      ! `scale` must be `TRUE` or `FALSE`, not a string.

---

    Code
      lwpls_fit(x, y, num_comps = 2)
    Condition
      Error in `lwpls_fit()`:
      ! `...` must be empty.
      x Problematic argument:
      * num_comps = 2

---

    Code
      fit <- lwpls_fit(x, y, neighbors = 500)
    Condition
      Warning:
      ! `neighbors` = 500 is larger than the number of training samples (90).
      i All 90 samples will be used.

# bad data are rejected

    Code
      lwpls_fit(list(1))
    Condition
      Error in `lwpls_fit()`:
      ! `lwpls_fit()` is not defined for a list.

---

    Code
      lwpls_fit(matrix(rnorm(20), 10), y[1:10])
    Condition
      Error in `lwpls_fit()`:
      ! The predictors must have unique, non-empty column names.

---

    Code
      lwpls_fit(iris[5:4], iris$Sepal.Length)
    Condition
      Error in `lwpls_fit()`:
      ! All predictors must be numeric.
      x Non-numeric predictor: `Species`.
      i Use the formula interface or a recipe to create dummy variables.

---

    Code
      lwpls_fit(x[1, ], y[1])
    Condition
      Error in `lwpls_fit()`:
      ! At least two training samples are required.

---

    Code
      lwpls_fit(x[0], y)
    Condition
      Error in `lwpls_fit()`:
      ! At least one predictor is required.

---

    Code
      lwpls_fit(x_na, y)
    Condition
      Error in `lwpls_fit()`:
      ! The predictors contain missing or infinite values.

---

    Code
      lwpls_fit(x, y_na)
    Condition
      Error in `lwpls_fit()`:
      ! The outcome contains missing or infinite values.

---

    Code
      lwpls_fit(x, factor(ifelse(y_na > 1, "a", "b")))
    Condition
      Error in `lwpls_fit()`:
      ! The outcome contains missing values.

---

    Code
      lwpls_fit(x, as.character(y > 1))
    Condition
      Error in `lwpls_fit()`:
      ! Character outcomes are not supported.
      i Convert the outcome to a factor for classification.

---

    Code
      lwpls_fit(x, factor(rep("a", nrow(x))))
    Condition
      Error in `lwpls_fit()`:
      ! The outcome factor must have at least two levels.

---

    Code
      lwpls_fit(x, data.frame(a = y, b = factor(y > 1)))
    Condition
      Error in `lwpls_fit()`:
      ! Unsupported outcome type.
      i Use a single factor (classification) or one or more numeric columns (regression).

# bad prediction arguments are rejected

    Code
      predict(fit, test[x_names], type = "class")
    Condition
      Error in `predict()`:
      ! `type` must be one of "numeric" or "raw", not "class".

---

    Code
      predict(fit, test[x_names], num_comp = 0)
    Condition
      Error in `predict()`:
      ! `num_comp` must be a single whole number >= 1.

---

    Code
      predict(fit, test[x_names], num_comp = 1:2)
    Condition
      Error in `predict()`:
      ! `num_comp` must be a single whole number >= 1.

---

    Code
      predict(fit, test[x_names], foo = 1)
    Condition
      Error in `predict()`:
      ! `...` must be empty.
      x Problematic argument:
      * foo = 1

