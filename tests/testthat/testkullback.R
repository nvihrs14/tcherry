context("kullback")
library(tcherry)

p_target <- array(c(0.1, 0.2, 0.05, 0.4, 0.005, 0.015, 0.13, 0.1),
                  dim = c(2, 2, 2),
                  dimnames = list("V1" = c("a", "b"),
                                  "V2" = c("a", "b"),
                                  "V3" = c("a", "b")))

p_approx <- array(c(0.01, 0.23, 0.06, 0.4, 0.005, 0.018, 0.137, 0.14),
                  dim = c(2, 2, 2),
                  dimnames = list("V1" = c("a", "b"),
                                  "V2" = c("a", "b"),
                                  "V3" = c("a", "b")))

KL <- 0.1 * log2(0.1 / 0.01) + 0.2 * log2(0.2 / 0.23) +
  0.05 * log2(0.05 / 0.06) + 0.4 * log2(0.4 / 0.4) +
  0.005 * log2(0.005 / 0.005) + 0.015 * log2(0.015 / 0.018) +
  0.13 * log2(0.13 / 0.137) + 0.1 * log2(0.1 / 0.14)

p_data_frame <- as.data.frame(p_target)

p_target_noname <- array(c(0.1, 0.2, 0.05, 0.4, 0.005, 0.015,
                           0.13, 0.1), dim = c(2, 2, 2))

p_target_names <- array(c(0.1, 0.2, 0.05, 0.4, 0.005,
                          0.015, 0.13, 0.1),
                        dim = c(2, 2, 2),
                        dimnames = list("V1" = NULL,
                                        "V2" = NULL,
                                        "V3" = NULL))

p_target_dimnames <- array(c(0.1, 0.2, 0.05, 0.4, 0.005, 0.015,
                             0.13, 0.1),
                           dim = c(2, 2, 2),
                           dimnames = list(c("a", "b"),
                                           c("a", "b"),
                                           c("a", "b")))

p_target_no_prop <- p_target
p_target_no_prop[1, 2, 2] <- 0.8

p_target_0 <- p_target
p_target_0[1, 2, 2] <- p_target_0[1, 2, 2] + p_target_0[2, 1, 1]
p_target_0[2, 1, 1] <- 0

p_approx_2 <- array(c(0.01, 0.23, 0.06, 0.4, 0.005, 0.018,
                      0.137, 0.14),
                    dim = c(2, 2, 2),
                    dimnames = list("A" = c("a", "b"),
                                    "B" = c("a", "b"),
                                    "C" = c("a", "b")))

p_approx_3 <- array(c(0.01, 0.23, 0.06, 0.4, 0.005, 0.018,
                      0.133, 0.14, 0.001, 0.001, 0.001, 0.001),
                    dim = c(2, 3, 2),
                    dimnames = list("V1" = c("a", "b"),
                                    "V2" = c("a", "b", "c"),
                                    "V3" = c("a", "b")))

test_that("kullback error messages work", {
  expect_error(kullback(p_data_frame, p_approx),
               "The input is not arrays.")
  expect_error(kullback(p_target_noname, p_approx),
               "The arrays must be named.")
  expect_error(kullback(p_target_dimnames, p_approx),
               "Variable names are missing from the arrays.")
  expect_error(kullback(p_target_names, p_approx),
               "Names of the variable levels are missing.")
  expect_error(kullback(p_target_no_prop, p_approx),
               "At least one of the given arrays is not a probability distribution.")
  expect_error(kullback(p_target_0, p_approx),
               "Some probabilities are zero.")
  expect_error(kullback(p_target, p_approx_2),
               "Distributions are not over the same universe.")
  expect_error(kullback(p_target, p_approx_3),
               "Distributions are not over the same universe.")
})

test_that("Kullback-Leibler is correct", {
  expect_equal(kullback(p_target, p_approx), KL)
  expect_equal(kullback(p_target, p_target), 0)
})
