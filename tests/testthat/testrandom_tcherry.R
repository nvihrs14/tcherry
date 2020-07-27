# https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17494
suppressWarnings(RNGversion("3.5.0"))

context("random_tcherry")
library(tcherry)

test_that("error messages work", {
  expect_error(random_tcherry("a", c(2, 3)), "n must be a single integer.")
  expect_error(random_tcherry(c(1, 2), c(2, 3)),
               "n must be a single integer.")
  expect_error(random_tcherry(2.3, c(2, 3)),
               "n must be a positive integer and at least 2.")
  expect_error(random_tcherry(1, c(2, 3)),
               "n must be a positive integer and at least 2.")
  expect_error(random_tcherry(- 1, c(2, 3)),
               "n must be a positive integer and at least 2.")
  expect_error(random_tcherry(2, c(2, 3), "a"),
               "noise must be a single numeric number.")
  expect_error(random_tcherry(2, c(2, 3), c(1, 1, 2)),
               "noise must be a single numeric number.")
  expect_error(random_tcherry(2, c(2, 3), - 3),
               "noise must be non-negative.")
  expect_error(random_tcherry(2, c("a", 3)),
               "n_levels must be a numeric vector.")
  expect_error(random_tcherry(2, matrix(1, 2, 2)),
               "n_levels must be a numeric vector.")
  expect_error(random_tcherry(2, c(2.1, 3)),
               "n_levels must be all positive integers.")
  expect_error(random_tcherry(2, c(- 2, 3)),
               "n_levels must be all positive integers.")
  expect_error(random_tcherry(2, c(0, 3)),
               "n_levels must be all positive integers.")
  expect_error(random_tcherry(2, c(2, 3, 3)),
               "The number of entries in n_level must be n.")
})

tch <- random_tcherry(5, rep(2, 5))

test_that("Check that the CPT sums to 1", {
  expect_true(sum(tch$CPTs[[1]]) == 1)
  expect_true(sum(tch$CPTs[[2]]) == 1)
  expect_true(all(colSums(tch$CPTs[[3]]) == 1))
  expect_true(all(colSums(tch$CPTs[[4]]) == 1))
  expect_true(all(colSums(tch$CPTs[[5]]) == 1))
})
