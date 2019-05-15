isDAG <- function(adj_matrix){
  
  while (any(rowSums(adj_matrix) == 0)) {
    idx <- which(rowSums(adj_matrix) == 0)[1]
    adj_matrix <- as.matrix(adj_matrix[- idx, - idx])
  }
  
  if (nrow(adj_matrix) == 0){
    res <- TRUE
  } else {
    res <- FALSE
    }
  res
}

#' Estimate conditional probability tables
#'
#' @description Estimates the conditional probability tables for
#' bayesian network models, where the structure is given by an
#' adjacency matrix.
#'
#' @param adj_matrix The adjacency matrix for the DAG.
#' @param data The data the probabilities should be estimated from.
#' @param bayes_smooth The additional cell counts for
#' bayesian estimation of probability tables.
#'
#' @return A list of the conditional probability tables for the
#' bayesian network. If the \code{bayes_smooth} argument is zero,
#' it is the maximum likelihood estimates. Otherwise, it is bayesian
#' estimates.
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 50, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 50, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 50, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 50, replace = TRUE))
#' 
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4))
#'
#' adj_matrix_DAG <- matrix(c(0, 0, 0, 0,
#'                            1, 0, 0, 0,
#'                            1, 0, 0, 0,
#'                            0, 1, 0, 0),
#'                           nrow = 4)
#'                           
#' rownames(adj_matrix_DAG) <- colnames(adj_matrix_DAG) <- names(data)
#'                           
#' CPT(adj_matrix_DAG, data)
#' CPT(adj_matrix_DAG, data, bayes_smooth = 1)
#' @export

CPT <- function(adj_matrix, data, bayes_smooth = 0){

  if (any(is.na(data))){
    warning(paste("The data contains NA values.",
                  "Theese will be excluded from tables,",
                  "which may be problematic.",
                  "It is highly recommended to manually take",
                  "care of NA values before using the data as input.",
                  sep = " "))
  }

  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }
  
  data <- data.frame(data, stringsAsFactors = FALSE)
  
  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }
  
  if (! is.matrix(adj_matrix)){
    stop("adj_matrix must be a matrix.")
  }
  
  if (any(diag(adj_matrix) == 1)){
    stop("The graph represented by adj_matrix contains loops.")
  }
  
  if (! is.numeric(adj_matrix)){
    stop("adj_matrix must be numeric.")
  }
  
  if (any(! c(adj_matrix) %in% 0:1)){
    stop(paste("adj_matrix must be an adjacency matrix for an unweighted graph.",
               "Therefore all entries must be 0 or 1.", sep = " "))
  }
  
  if (is.null(colnames(adj_matrix)) | is.null(rownames(adj_matrix))){
    stop("adj_matrix must be named.")
  }
  
  if (any(colnames(adj_matrix) != rownames(adj_matrix))){
    stop("Names of columns and rows in adj_matrix must be the same.")
  }
  
  if (length(setdiff(colnames(adj_matrix), names(data))) != 0){
    stop("The names of adj_matrix must be variable names in data.")
  }
  
  if (! isDAG(adj_matrix)){
    stop("adj_matrix is not a DAG.")
  }
  
  if (length(bayes_smooth) > 1){
    stop("bayes_smooth must be a single non-negative value.")
  }
  else if (!is.numeric(bayes_smooth)) {
    stop("bayes_smooth must be numeric.")
  }
  else if (bayes_smooth < 0){
    stop("bayes_smooth must be a non-negative numeric value.")
  }
  
  nodes <- rownames(adj_matrix)
  FUN <- function(node){
    parents_idx <- which(adj_matrix[, node] == 1)
    parents <- nodes[parents_idx]
    tab <- table(data[, c(node, parents)]) + bayes_smooth
    
    if (length(parents) == 0){
      mar <- NULL
      names(dimnames(tab)) <- node
    } else {
      tab_parents <- table(data[, c(parents)]) + bayes_smooth
      if (any(tab_parents == 0)){
        stop("Some cell counts of parent configurations are zero.
             Consider using the bayes_smooth argument.")
      }
      mar <- (1:length(parents)) + 1
      }
    prop.table(tab, margin = mar)
  }

  CPT_list <- lapply(nodes, FUN)
  names(CPT_list) <- nodes
  CPT_list
}
