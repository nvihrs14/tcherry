% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CL.R
\name{ChowLiu}
\alias{ChowLiu}
\title{Determines the Chow-Liu tree for data}
\usage{
ChowLiu(data, root = NULL, bayes_smooth = 0, CPTs = TRUE, ...)
}
\arguments{
\item{data}{The data the tree structure should be based on.}

\item{root}{An optional argument, choosing a userspecified
root for the tree.}

\item{bayes_smooth}{Additional cell counts for bayesian
estimation of probability tables.}

\item{CPTs}{A logical value indicating whether conditional probability
tables should be estimated from data or not.}

\item{...}{Additional parameters passed to \code{MI2}.}
}
\value{
A list containing the following components:
\itemize{
\item \code{skeleton_adj} The adjacency matrix for the skeleton
of the Chow-Liu tree.
\item \code{adj_DAG} The adjacency matrix of the resulting DAG.
\item \code{CPTs} The estimated conditional probability tables
of the bayesian network if estimated. Otherwise the logical value
\code{FALSE}.
}
}
\description{
Determines the structure and the conditional
probability tables for the Chow-Liu tree fitted to data.
}
\examples{
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

CL <- ChowLiu(data, root = 'var1', smooth = 0.1)
}
\seealso{
\code{\link{MI2}} for mutual information,
\code{\link{CPT}} for probability tables and
\code{\link{is_acyclic}} for a test for cycles.
}
\author{
Katrine Kirkeby, \email{enir_tak@hotmail.com}

Maria Knudsen, \email{mariaknudsen@hotmail.dk}

Ninna Vihrs, \email{ninnavihrs@hotmail.dk}
}
