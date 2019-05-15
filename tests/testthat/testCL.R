context("is_acyclic")
library(tcherry)

adj_matrix_cyclic <- matrix(c(0, 1, 1, 1,
                              1, 0, 0, 1,
                              1, 0, 0, 0,
                              1, 1, 0, 0),
                            nrow = 4)

adj_matrix_acyclic <- matrix(c(0, 0, 1, 1,
                               0, 0, 0, 1,
                               1, 0, 0, 0,
                               1, 1, 0, 0),
                             nrow = 4)

adj_matrix_loop <- matrix(c(1, 0, 1, 1,
                            0, 0, 0, 1,
                            1, 0, 0, 0,
                            1, 1, 0, 0),
                          nrow = 4)

adj_matrix_acyclic_disconnected <- matrix(c(0, 0, 1, 1,
                                            0, 0, 0, 0,
                                            1, 0, 0, 0,
                                            1, 0, 0, 0),
                                          nrow = 4)

adj_no_num <- matrix(c("0", 0, 1, 1,
                        0, 0, 0, 0,
                        1, 0, 0, 0,
                        1, 0, 0, 0),
                      nrow = 4)

adj_non_adj <- matrix(c(0, 0, 2, 1,
                        0, 0, 0, 0,
                        2, 0, 0, 0,
                        1, 0, 0, 0),
                      nrow = 4)

adj_no_sym <- matrix(c(0, 0, 1, 1,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       1, 0, 0, 0),
                      nrow = 4)

test_that("is_acyclic is working", {
  expect_error(is_acyclic(1), "Argument must be a matrix.")
  expect_error(is_acyclic(adj_matrix_loop),
               "The graph represented by the matrix contains loops.")
  expect_error(is_acyclic(adj_no_num), "Argument must be numeric.")
  expect_error(is_acyclic(adj_non_adj), 
               paste("Argument must be an adjacency matrix for an unweighted graph.",
                     "Therefore all entries must be 0 or 1.", sep = " "))
  expect_error(is_acyclic(adj_no_sym), 
               paste("Only undirected graphs are supported so argument must be",
                     "symmetric. This includes that rownames must equal colnames.", 
                     sep = " "))
  expect_true(is_acyclic(adj_matrix_acyclic))
  expect_true(is_acyclic(adj_matrix_acyclic_disconnected))
  expect_false(is_acyclic(adj_matrix_cyclic))
})

context("CPT")

set.seed(43)
var1 <- c(sample(c(1, 2), 50, replace = TRUE))
var2 <- var1 + c(sample(c(1, 2), 50, replace = TRUE))
var3 <- var1 + c(sample(c(0, 1), 50, replace = TRUE,
                        prob = c(0.9, 0.1)))
var4 <- c(sample(c(1, 2), 50, replace = TRUE))

data <- data.frame("var1" = as.character(var1),
                   "var2" = as.character(var2),
                   "var3" = as.character(var3),
                   "var4" = as.character(var4))

data_num <- data.frame("var1" = as.character(var1),
                       "var2" = as.character(var2),
                       "var3" = as.character(var3),
                       "var4" = var4)

adj_matrix_DAG <- matrix(c(0, 0, 0, 0,
                           1, 0, 0, 0,
                           1, 0, 0, 0,
                           0, 1, 0, 0),
                           nrow = 4)

rownames(adj_matrix_DAG) <- colnames(adj_matrix_DAG) <- names(data)

adj_matrix_loop <- matrix(c(1, 0, 0, 0,
                            1, 0, 0, 0,
                            1, 0, 0, 0,
                            0, 1, 0, 0),
                          nrow = 4)

adj_matrix_no_num <- matrix(c(0, "0", 0, 0,
                              1, 0, 0, 0,
                              1, 0, 0, 0,
                              0, 1, 0, 0),
                            nrow = 4)

adj_matrix_no_adj <- matrix(c(0, 2, 0, 0,
                              1, 0, 0, 0,
                              1, 0, 0, 0,
                              0, 1, 0, 0),
                            nrow = 4)

adj_matrix_no_name <- matrix(c(0, 0, 0, 0,
                               1, 0, 0, 0,
                               1, 0, 0, 0,
                               0, 1, 0, 0),
                             nrow = 4)

adj_matrix_wrong_name <- matrix(c(0, 0, 0, 0,
                                  1, 0, 0, 0,
                                  1, 0, 0, 0,
                                  0, 1, 0, 0),
                                nrow = 4)

rownames(adj_matrix_wrong_name) <- names(data)
colnames(adj_matrix_wrong_name) <- rownames(adj_matrix_wrong_name)[4:1]

adj_matrix_wrong_name2 <- matrix(c(0, 0, 0, 0,
                                   1, 0, 0, 0,
                                   1, 0, 0, 0,
                                   0, 1, 0, 0),
                                 nrow = 4)

rownames(adj_matrix_wrong_name2) <- colnames(adj_matrix_wrong_name2) <- letters[1:4]

adj_matrix_no_DAG <- matrix(c(0, 0, 0, 1,
                              1, 0, 0, 0,
                              1, 0, 0, 0,
                              0, 1, 0, 0),
                            nrow = 4)

rownames(adj_matrix_no_DAG) <- colnames(adj_matrix_no_DAG) <- names(data)

CPT_var2_var1 <- c(0.6, 0.4, 0, 0, 0.35, 0.65)
CPT_ex <- CPT(adj_matrix_DAG, data)

test_that("CPT's are calculated correctly", {
  expect_equal(c(CPT_ex$var2), CPT_var2_var1)
  expect_equal(sum(CPT_ex$var3[, 1]), 1)
  expect_equal(sum(CPT_ex$var3[, 2]), 1)
  expect_equal(sum(CPT_ex$var4[, 1]), 1)
  expect_equal(sum(CPT_ex$var1), 1)
})

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(CPT(adj_matrix_DAG, data_na),
                 paste("The data contains NA values.",
                       "Theese will be excluded from tables,",
                       "which may be problematic.",
                       "It is highly recommended to manually take",
                       "care of NA values before using the data as input.",
                       sep = " "))
})

test_that("Error messages work", {
  expect_error(CPT(adj_matrix_DAG, c(1, 2)), 
               "data must be a data frame or a matrix.")
  expect_error(CPT(adj_matrix_DAG, data_num),
               "Some columns are not characters or factors.")
  expect_error(CPT("var1", data),
               "adj_matrix must be a matrix.")
  expect_error(CPT(adj_matrix_loop, data),
               "The graph represented by adj_matrix contains loops.")
  expect_error(CPT(adj_matrix_no_num, data), "adj_matrix must be numeric.")
  expect_error(CPT(adj_matrix_no_adj, data),
            paste("adj_matrix must be an adjacency matrix for an unweighted graph.",
                     "Therefore all entries must be 0 or 1.", sep = " "))
  expect_error(CPT(adj_matrix_no_name, data),
               "adj_matrix must be named.")
  expect_error(CPT(adj_matrix_wrong_name, data),
               "Names of columns and rows in adj_matrix must be the same.")
  expect_error(CPT(adj_matrix_wrong_name2, data),
               "The names of adj_matrix must be variable names in data.")
  expect_error(CPT(adj_matrix_no_DAG, data),"adj_matrix is not a DAG.")
  expect_error(CPT(adj_matrix_DAG, data, bayes_smooth = c(1, 2)),
               "bayes_smooth must be a single non-negative value.")
  expect_error(CPT(adj_matrix_DAG, data, bayes_smooth = "a"), 
               "bayes_smooth must be numeric.")
  expect_error(CPT(adj_matrix_DAG, data, bayes_smooth = - 1),
               "bayes_smooth must be a non-negative numeric value.")
  })

context("ChowLiu")

data_matrix <- as.matrix(data)

data_numeric <- data.frame("var1" = as.character(var1),
                           "var2" = as.character(var2),
                           "var3" = as.character(var3),
                           "var4" = var4)

test_that("Input is specified correctly", {
  expect_error(ChowLiu(data, root = "var5"),
               "The specified root is not a node.")
  expect_error(ChowLiu(data_numeric),
               "Some columns are not characters or factors.")
  expect_error(ChowLiu(var1),
               "data must be a data frame or a matrix.")
  expect_error(ChowLiu(data, bayes_smooth = c(1, 2)),
               "bayes_smooth must be a single non-negative value.")
  expect_error(ChowLiu(data, bayes_smooth = "a"), "bayes_smooth must be numeric.")
  expect_error(ChowLiu(data, bayes_smooth = - 1),
               "bayes_smooth must be a non-negative numeric value.")
})

adj_matrix_skeleton <- matrix(c(0, 1, 1, 0,
                                1, 0, 0, 1,
                                1, 0, 0, 0,
                                0, 1, 0, 0),
                              nrow = 4)

rownames(adj_matrix_skeleton) <-
  colnames(adj_matrix_skeleton) <- names(data)

CL <- ChowLiu(data, root = "var1", smooth = 0.1)
CL_mat <- ChowLiu(data_matrix, root = "var1", smooth = 0.1)

test_that("Adjacency matrices are correct", {
  expect_equal(CL$adj_DAG, adj_matrix_DAG)
  expect_equal(CL$skeleton_adj, adj_matrix_skeleton)
  expect_equal(CL_mat$adj_DAG, adj_matrix_DAG)
})

data_na <- data
data_na[1, 1] <- NA

test_that("Warning message works", {
  expect_warning(ChowLiu(data_na, smooth = 0.1))
})