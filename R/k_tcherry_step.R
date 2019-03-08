#' Determine a k'th order t-cherry tree from data
#'
#' @description Determine the structure of a k'th order t-cherry tree
#' from data based on a greedy stepwise approach.
#'
#' @param data The data the tree structure should be based on.
#' @param k The order of the t-cherry tree.
#' @param ... Additional arguments passed to \code{MIk}.
#'
#' @details Notice that for \eqn{k=3} it is the same as using
#' \code{tcherry_step} and for \eqn{k=2} it is the same as using
#' \code{ChowLiu}.
#'
#' The algorithm for constructing the t-cherry tree from
#' data is based on an atempt to minimize the Kullback-Leibler
#' divergence. The first cherry is chosen as the k variables with
#' highest mutual information. This is the preliminary t-cherry
#' tree. Then all possible new cherries are added stepwise to this
#' tree and the weight \deqn{\sum MI(Clique) - \sum MI(Separator)} is
#' calculated.
#' The first sum is over the cliques and the second over the
#' separators of the junction tree of the preliminary t-cherry tree.
#' The one with the highest weight is chosen as the new preliminary
#' t-cherry tree, and the procedure is repeated untill all variables
#' has been added.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{adj_matrix} The adjacency matrix for the k'th order
#' t-cherry tree.
#' \item \code{weight} The weight of the final k'th order t-cherry tree.
#' \item \code{cliques} A list containing the cliques (cherries) of
#'  the k'th order t-cherry tree.
#' \item \code{separators} A list containing the separators of a
#' junction tree for the k'th order t-cherry tree.
#' }
#'
#' @author
#' Katrine Kirkeby, \email{enir_tak@@hotmail.com}
#'
#' Maria Knudsen, \email{mariaknudsen@@hotmail.dk}
#'
#' Ninna Vihrs, \email{ninnavihrs@@hotmail.dk}
#'
#' @seealso \code{\link{MIk}} for mutual
#' information of k variables.
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
#' # smooth used in MIk
#' (tch <- k_tcherry_step(data, 3, smooth = 0.1))
#'
#' # For plotting
#' library(gRbase)
#' library(Rgraphviz)
#' tcherry_tree <- as(tch$adj_matrix, "graphNEL")
#' plot(tcherry_tree)
#'
#' # For probability propagation
#' library(gRain)
#' model <- grain(tcherry_tree, data = data, smooth = 0.1)
#' querygrain(model)
#' @export

k_tcherry_step <- function(data, k, ...){
  if (! (is.data.frame(data) | is.matrix(data))) {
    stop("data must be a data frame or a matrix.")
  }

  if (! all(sapply(data, function(x){
    is.character(x) | is.factor(x)
  }
  ))){
    stop("Some columns are not characters or factors.")
  }

  if (length(k) != 1){
    stop("k must be a single positive integer.")
  }

  if (k %% 1 != 0 | k <= 1){
    stop("k must be a positive integer and at least 2.")
  }

  nodes <- colnames(data)
  n_var <- length(nodes)
  adj_matrix <- matrix(0, nrow = n_var, ncol = n_var,
                       dimnames = list(nodes, nodes))
  cliques <- as.list(rep(NA, n_var - (k - 1)))
  separators <- as.list(rep(NA, n_var - k))

  # Adding first cherry
  poss_cliq <- utils::combn(nodes, k)
  poss_cliq <- split(poss_cliq, rep(1:ncol(poss_cliq),
                                    each = nrow(poss_cliq)))

  MI <- sapply(poss_cliq, MIk, data = data, ...)
  idx.max <- which.max(MI)

  tcherry_nodes <- poss_cliq[[idx.max]]
  cliques[[1]] <- tcherry_nodes
  nodes_remaining <- setdiff(nodes, tcherry_nodes)
  tcherry_hyperedges <- utils::combn(tcherry_nodes, k-1)
  tcherry_hyperedges <- split(tcherry_hyperedges,
                              rep(1:ncol(tcherry_hyperedges),
                                  each = nrow(tcherry_hyperedges)))

  adj_matrix[tcherry_nodes, tcherry_nodes] <- 1
  diag(adj_matrix[tcherry_nodes, tcherry_nodes]) <- 0

  weight <- max(MI)

  if (length(nodes_remaining) != n_var){

  idx_list <- 1

  n_nodes_remaining_median <- floor((n_var - k) / 2 + 1)
  n_hyp_edges_median <- 1 + (k - 1) * ((n_var - n_nodes_remaining_median)
                                       - (k - 1))

  weight_cliq_sep <- MI_cliq <- MI_sep <- new_var <-
    rep(NA, n_nodes_remaining_median * n_hyp_edges_median)

  new_cliques_list <- new_seps_list <-
    as.list(weight_cliq_sep)

  dat_new_poss <- data.frame(new_cliq = I(new_cliques_list),
                             new_sep = I(new_seps_list),
                             new_var = new_var,
                             MI_cliq = MI_cliq,
                             MI_sep = MI_sep,
                             weight_increase = weight_cliq_sep)

  idx.dat <- 1
  for (i in 1:length(tcherry_hyperedges)){
    for (var in nodes_remaining){
      dat_new_poss$new_sep[[idx.dat]] <- tcherry_hyperedges[[i]]
      dat_new_poss$new_cliq[[idx.dat]] <-
        c(dat_new_poss$new_sep[[idx.dat]], var)

      dat_new_poss$MI_sep[idx.dat] <-
        MIk(dat_new_poss$new_sep[[idx.dat]], data, ...)
      dat_new_poss$MI_cliq[idx.dat] <-
        MIk(dat_new_poss$new_cliq[[idx.dat]], data, ...)

      dat_new_poss$weight_increase[idx.dat] <-
        dat_new_poss$MI_cliq[idx.dat] - dat_new_poss$MI_sep[idx.dat]

      dat_new_poss$new_var[[idx.dat]] <- var
      idx.dat <- idx.dat + 1
    }
  }

  idx_max_weight <- which.max(dat_new_poss$weight_increase)
  weight <- weight + dat_new_poss$weight_increase[idx_max_weight]

  new_clique <- dat_new_poss$new_cliq[[idx_max_weight]]
  new_sep <- dat_new_poss$new_sep[[idx_max_weight]]
  new_var <- dat_new_poss$new_var[idx_max_weight]

  tcherry_nodes <- c(tcherry_nodes, new_var)
  nodes_remaining <- setdiff(nodes, tcherry_nodes)

  cliques[[2]] <- new_clique
  separators[[1]] <- new_sep
  idx.list <- 2

  new_hyper_edges <- utils::combn(new_clique, k - 1)
  new_hyper_edges <- split(new_hyper_edges,
                           rep(1:ncol(new_hyper_edges),
                               each = nrow(new_hyper_edges)))
  idx.new <- sapply(new_hyper_edges, function(e){
    length(setdiff(e, new_var)) != k - 1
  })
  new_hyper_edges <- new_hyper_edges[idx.new]

  adj_matrix[new_clique, new_clique] <- 1
  diag(adj_matrix[new_clique, new_clique]) <- 0

  idx.new.var <- dat_new_poss$new_var == new_var &
    !is.na(dat_new_poss$new_var)
  dat_new_poss[idx.new.var, ] <- NA
  ord <- order(dat_new_poss$weight_increase, decreasing = TRUE)
  dat_new_poss <- dat_new_poss[ord, ]
  }

  while (length(tcherry_nodes) != n_var) {
    idx.dat <- which(is.na(dat_new_poss$weight_increase))[1]
    for (i in 1:length(new_hyper_edges)){
      for (var in nodes_remaining){
        dat_new_poss$new_sep[[idx.dat]] <- new_hyper_edges[[i]]
        dat_new_poss$new_cliq[[idx.dat]] <-
          c(dat_new_poss$new_sep[[idx.dat]], var)

        dat_new_poss$MI_sep[idx.dat] <-
          MIk(dat_new_poss$new_sep[[idx.dat]], data, ...)
        dat_new_poss$MI_cliq[idx.dat] <-
          MIk(dat_new_poss$new_cliq[[idx.dat]], data, ...)

        dat_new_poss$weight_increase[idx.dat] <-
          dat_new_poss$MI_cliq[idx.dat] - dat_new_poss$MI_sep[idx.dat]

        dat_new_poss$new_var[[idx.dat]] <- var
        idx.dat <- idx.dat + 1
      }
    }

    idx_max_weight <- which.max(dat_new_poss$weight_increase)
    weight <- weight + dat_new_poss$weight_increase[idx_max_weight]

    new_clique <- dat_new_poss$new_cliq[[idx_max_weight]]
    new_sep <- dat_new_poss$new_sep[[idx_max_weight]]
    new_var <- dat_new_poss$new_var[idx_max_weight]

    tcherry_nodes <- c(tcherry_nodes, new_var)
    nodes_remaining <- setdiff(nodes, tcherry_nodes)

    cliques[[idx.list + 1]] <- new_clique
    separators[[idx.list]] <- new_sep
    idx.list <- idx.list + 1

    new_hyper_edges <- utils::combn(new_clique, k - 1)
    new_hyper_edges <- split(new_hyper_edges,
                             rep(1:ncol(new_hyper_edges),
                                 each = nrow(new_hyper_edges)))
    idx.new <- sapply(new_hyper_edges, function(e){
      length(setdiff(e, new_var)) != k - 1
    })
    new_hyper_edges <- new_hyper_edges[idx.new]

    adj_matrix[new_clique, new_clique] <- 1
    diag(adj_matrix[new_clique, new_clique]) <- 0

    idx.new.var <- dat_new_poss$new_var == new_var &
      !is.na(dat_new_poss$new_var)
    dat_new_poss[idx.new.var, ] <- NA
    ord <- order(dat_new_poss$weight_increase, decreasing = TRUE)
    dat_new_poss <- dat_new_poss[ord, ]


  }

  return(list("adj_matrix" = adj_matrix,
              "weight" = weight,
              "cliques" = cliques,
              "separators" = separators))
}
