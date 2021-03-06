#' Determine a (k + 1)'th order t-cherry tree from a k'th order t-cherry
#' tree
#'
#' @description Determine the structure of a (k + 1)'th order
#' t-cherry tree from a k'th order t-cherry tree.
#'
#' @param tch_cliq A list containing the cliques of the k'th order
#' t-cherry tree.
#' @param data The data the tree structure should be based on.
#' @param ... Additional arguments passed to \code{MIk}.
#'
#' @details The algorithms for constructing the (k + 1)'th order t-cherry
#' tree from a k'th order t-cherry tree are greedy algorithms.
#' \code{increase_order1} attempts to maximize the sum of
#' mutual information of the cliques and
#' \code{increase_order2} attempts to maximize the weight of
#' the junction tree. \code{increase_order2} is also a faster
#' implementation and a faster alternative to creating a (k + 1)'th order
#' t-cherry tree directly from data with \code{k_tcherry_step}. It is
#' therefore recommended to use this one, and \code{increase_order1}
#' is primarily kept for historical reasons.
#'
#' In \code{increase_order1} the procedure is:
#' \itemize{
#' \item Starting from the k'th order t-cherry tree, make a complete set
#' of the (k + 1) variables with highest mutual information which
#' satisfies
#' that this only adds one edge to the original graph. Remove these
#' (k + 1)
#' variables from later consideration.
#' \item Make a complete set of the (k + 1) variables with highest mutual
#' information which satisfies that this only adds one edge to the graph
#' and that k of the variables are in earlier created cliques of size
#' (k + 1). Remove these (k + 1) variables from later consideration.
#' \item Continue until all variables are in a clique of size (k + 1).
#' }
#'
#' For \code{increase_order2} the procedure is: Start with
#' the k'th order t-cherry tree T_k and set T_(k + 1) = T_k.
#' Choose the set of (k + 1) variables with highest mutual information which
#' satisfies that making the set complete in T_(k + 1) adds only one edge.
#' Let this set be the first cherry/clique and make it complete in T_(k + 1).
#' Consider all possible new cherries of size (k + 1).
#' A new cherry is possible if k of the variables are already in an
#' existing cherry, and making the set complete only adds one edge in T_(k + 1).
#' For each new cherry, calculate the weight \deqn{MI(clique) - MI(separator)}
#' for the new clique and separator of the junction tree for the preliminary
#' (k + 1)'th order t-cherry tree. Add the cherry with the highest weight to
#' T_(k + 1). Repeat the procedure until T_(k + 1) is a (k + 1)'th order
#' t-cherry tree. 
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_matrix} The adjacency matrix for the (k + 1)'th order
#' t-cherry tree.
#' \item \code{weight} Weight of the junction tree (only for \code{increase_order2})
#' \item \code{cliques} A list containing the cliques (cherries) of
#'  the (k + 1)'th order t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the (k + 1)'th order t-cherry tree.
#' \item \code{n_edges} The number of edges in the resulting graph.
#' }
#' 
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{MIk}} for mutual information of k variables.
#'
#' @examples
#' set.seed(43)
#' var1 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var2 <- var1 + c(sample(c(1, 2), 100, replace = TRUE))
#' var3 <- var1 + c(sample(c(0, 1), 100, replace = TRUE,
#'                         prob = c(0.9, 0.1)))
#' var4 <- c(sample(c(1, 2), 100, replace = TRUE))
#' var5 <- var2 + var3
#' var6 <- var1 - var4 + c(sample(c(1, 2), 100, replace = TRUE))
#' var7 <- c(sample(c(1, 2), 100, replace = TRUE))
#'
#' data <- data.frame("var1" = as.character(var1),
#'                    "var2" = as.character(var2),
#'                    "var3" = as.character(var3),
#'                    "var4" = as.character(var4),
#'                    "var5" = as.character(var5),
#'                    "var6" = as.character(var6),
#'                    "var7" = as.character(var7))
#'
#' ChowLiu_cliques <- list(c("var1", "var5"),
#'                         c("var2", "var5"),
#'                         c("var3", "var5"),
#'                         c("var3", "var7"),
#'                         c("var4", "var6"),
#'                         c("var5", "var6"))
#'                         
#' # smooth used in MIk
#' (tch <- increase_order1(ChowLiu_cliques, data, smooth = 0.1))
#' (tch2 <- increase_order2(ChowLiu_cliques, data, smooth = 0.1))
#' @export

increase_order1 <- function(tch_cliq, data, ...){

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

  data <- as.data.frame(data)

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (! is.list(tch_cliq)){
    stop("Cliques must be given in a list, each entry containing a vector
         with the names of the variables in the clique.")
  }

  if (! compare::compare(unique(unlist(tch_cliq)), colnames(data),
                       ignoreOrder = TRUE)$result){
    stop("The column names of data must be the same as the variable
         names in tch_cliq. All variables in data must be in at least
         one clique.")
  }

  if (length(unique(sapply(tch_cliq, length))) != 1){
    stop("tch_cliq should be the cliques of a k'th order t-cherry tree.
         Therefore they should all have the same length k.")
  }

  nodes <- colnames(data)
  n_var <- length(nodes)
  n_cliq <- length(tch_cliq)
  k <- length(tch_cliq[[1]])

  if (n_var < (k + 1)){
    stop("It takes at least k plus 1 variables to fit a k plus 1'th order t-cherry tree.")
  }

  tch_adj <- matrix(0, nrow = n_var, ncol = n_var)
  rownames(tch_adj) <- colnames(tch_adj) <- nodes
  for (i in 1:n_cliq) {
    tch_adj[tch_cliq[[i]], tch_cliq[[i]]] <- 1
    diag(tch_adj[tch_cliq[[i]], tch_cliq[[i]]]) <- 0
  }

  if (! all(gRbase::triangulateMAT(tch_adj) == tch_adj)){
    stop("The cliques do not come from a triangulated graph.
         The cliques should correspond to a k'th order t-cherry tree
         so it must be triangulated.")
  }

  if (sum(tch_adj) / 2 != (k - 1) * n_var - (1 / 2) * (k - 1) * k){
    stop("The graph corresponding to the cliques does not have the
         correct number of edges for a k'th order t-cherry tree.")
  }

  k_plus_1_pairs <- utils::combn(nodes, k + 1)
  k_plus_1_pairs <- split(k_plus_1_pairs, rep(1:ncol(k_plus_1_pairs),
                                         each = nrow(k_plus_1_pairs)))

  MI <- sapply(k_plus_1_pairs, MIk, data = data, ...)

  ord_idx <- order(MI, decreasing = TRUE)
  MI <- MI[ord_idx]
  k_plus_1_pairs <- k_plus_1_pairs[ord_idx]

  n_edges <- sum(tch_adj) / 2
  tcherry_nodes <- c()

  # Adding the first cherry.
  i <- 1
  while (length(tcherry_nodes) == 0) {

    adj_matrix_temp <- tch_adj
    adj_matrix_temp[k_plus_1_pairs[[i]], k_plus_1_pairs[[i]]] <- 1
    diag(adj_matrix_temp[k_plus_1_pairs[[i]], k_plus_1_pairs[[i]]]) <- 0

    if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
      adj_matrix <- adj_matrix_temp
      tcherry_nodes <- k_plus_1_pairs[[i]]
      n_edges <- n_edges + 1
      MI <- MI[- i]
      k_plus_1_pairs <- k_plus_1_pairs[- i]
      cliques <- list(tcherry_nodes)
    }
    i <- i + 1
  }

  # Adding remaining nodes via new cherries.
  i <- l <- 1
  j <- 2

  separators <- list()

  while (length(tcherry_nodes) < n_var) {
    n_var_new <- length(setdiff(k_plus_1_pairs[[i]], tcherry_nodes))

    if (n_var_new == 1){

      adj_matrix_temp <- adj_matrix
      adj_matrix_temp[k_plus_1_pairs[[i]], k_plus_1_pairs[[i]]] <- 1
      diag(adj_matrix_temp[k_plus_1_pairs[[i]],
                           k_plus_1_pairs[[i]]]) <- 0

      if ((sum(adj_matrix_temp) / 2 ) == (n_edges + 1)){
        adj_matrix <- adj_matrix_temp
        cliques[[j]] <- k_plus_1_pairs[[i]]
        separators[[l]] <- intersect(cliques[[j]], tcherry_nodes)

        tcherry_nodes <- unique(c(tcherry_nodes, cliques[[j]]))
        n_edges <- n_edges + 1

        MI <- MI[- i]
        k_plus_1_pairs <- k_plus_1_pairs[- i]

        i <- 0
        j <- j + 1
        l <- l + 1
      }
    }

    if(n_var_new == 0){
      MI <- MI[- i]
      k_plus_1_pairs <- k_plus_1_pairs[- i]
      i <- i - 1
    }
    i <- i + 1
  }

  n_edges_graph <- sum(adj_matrix) / 2

  return(list("adj_matrix" = adj_matrix,
              "cliques" = cliques,
              "separators" = separators,
              "n_edges" = n_edges_graph))
}
