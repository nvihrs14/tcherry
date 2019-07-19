# tcherry

The package is meant for learning the structure of the type of graphical models called t-cherry trees from data. For more information about t-cherry trees see for instance the vignette or Kovács & Szántai (2012). The functions can only be used for categorical data. The purpose of the package is mainly to provide functions for learning the structure of the graph from given data with no missing values. The functions are attempting to find structures of maximal likelihood. If the corresponding graphical model is to be used in connection with probability propagation or prediction other packages are needed (see the vignette for examples). 

## Install

Without vignettes
`remotes::install_github("nvihrs14/tcherry")`

With vignettes
`remotes::install_github("nvihrs14/tcherry",  build_opts = c("--no-resave-data", "--no-manual"))`

If there are problems with viewing documentation or vignettes, it is recommended to restart the R session.

Note that the package requires the following R-packages, which are automatically installed with the package:
    Rdpack, utils, gRbase, compare, rGrain, Rgraphviz and stats.
    
The package gRbase further requires the packages graph and RBGL which may have to be installed from Bioconductor for instance with

`BiocManager::install(c("graph", "RBGL"))`

## Main functions (see vignette for more details)

-__`increase_order2`__: Determine a (k + 1)'th order t-cherry tree from a k'th order t-cherry tree by a greedy search.

-__`increase_order_complete_search`__: Determine the (k + 1)'th order t-cherry tree from a k'th order t-cherry tree by a complete search.

-__`k_tcherry_p_lookahead`__: Determine a k'th order t-cherry tree from data by adding p cliques at a time by a greedy search. Note that if p is the total number of cliques in a k'th order t-cherry tree with the desired number of vertices, this is a complete search.

-__`thin_edges`__: Thinning of edges in an undirected graphical model with a triangulated graph.

-__`compute_BIC_junction_tree`__: Calculates the BIC value for a graphical model from a junction tree for the graph.

## Example usage

To demonstrate the main functions in this package consider the car evaluation data set from UCI Machine Learning Repository (Dau & Graff 2017). This data set contains 7 variables (all categorical) with 1728 observations for each and no missing values. The variables are describing different aspects of the car such as the estimated safety of the car, the number of doors etc. To find a graphical structure of a third order t-cherry tree for this data the function `k_tcherry_p_lookahead` is used. It is chosen to add just one clique at a time in the greedy search procedure.

``` r
library(tcherry)
car <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data",
          header = FALSE, sep = ",", dec = ".")
names(car) <- c("buying", "maint", "doors", "persons", "lug_boot",
                  "safety", "class")
tch3 <- k_tcherry_p_lookahead(data = car, k = 3, p = 1, smooth = 0.001)
tch3$adj_matrix
#>            buying maint doors persons lug_boot safety class
#> buying        0     1     0       0        0      1     1
#> maint         1     0     0       0        0      0     1
#> doors         0     0     0       0        1      0     1
#> persons       0     0     0       0        0      1     1
#> lug_boot      0     0     1       0        0      1     1
#> safety        1     0     0       1        1      0     1
#> class         1     1     1       1        1      1     0

```

Note that the smooth argument is added to cell counts when estimating probabilities to avoid zero probabilities, which would make some calculations invalid. 

The graphical structure of af fourth order t-cherry tree for this data can be found by using the same function as above with k = 4. However in this case, it is chosen to show how `increase_order2` can be used to increase the order of the fitted third order t-cherry tree. The typical reason for this choice will be to save time, but often at the cost of a fitted structure of smaller likelihood.

``` r
tch4 <- increase_order2(tch_cliq = tch3$cliques, data = car, smooth = 0.001)
tch4$adj_matrix
#>            buying maint doors persons lug_boot safety class
#> buying        0     1     0       1        1      1     1
#> maint         1     0     0       0        0      1     1
#> doors         0     0     0       0        1      1     1
#> persons       1     0     0       0        0      1     1
#> lug_boot      1     0     1       0        0      1     1
#> safety        1     1     1       1        1      0     1
#> class         1     1     1       1        1      1     0

```

Note that the smooth argument is added for the same reasons as above, and the given third order t-cherry tree is represented by its cliques. 

It can now be attempted to simplify this model by thinning the edges.

``` r
tch_thinning <- thin_edges(cliques = tch4$cliques, separators = tch4$separators, data = car,
smooth = 0.001)
tch_thinning$adj_matrix
#>            buying class doors lug_boot maint persons safety
#> buying        0     1     0        0     1       0      1
#> class         1     0     0        1     1       1      1
#> doors         0     0     0        0     0       0      0
#> lug_boot      0     1     0        0     0       0      1
#> maint         1     1     0        0     0       0      0
#> persons       0     1     0        0     0       0      1
#> safety        1     1     0        1     0       1      0

tch_thinning$n_edges_removed
#> 6

```
Notice that in this function the structure is represented by the cliques and separators of its junction tree and that the resulting graph is no longer a t-cherry tree. However, the function will always return a triangulated graph. In this case, six edges has been deleted from the graph. The order in which the edges are deleted is random, so the resulting structure may vary.   

The three fitted structures can be compared by calculating a BIC score.

``` r
compute_BIC_junction_tree(cliques = tch3$cliques, separators = tch3$separators, data = car, smooth = 0.001)
#> -20079.89

compute_BIC_junction_tree(cliques = tch4$cliques, separators = tch4$separators, data = car, smooth = 0.001)
#> -21572.4

compute_BIC_junction_tree(cliques = tch_thinning$cliques, separators = tch_thinning$separators, data = car,
smooth = 0.001)
#> -19923.95
```
Because of the way this score is calculated, it is desired to get a high BIC score. The structure with the highest score is the thinned fourth order t-cherry tree. 

## For more help

See documentation included in package (vignettes and man) at <https://github.com/nvihrs14/tcherry>

## Feedback

If you experience any problems with this package, you are welcome to contact the authors by email. Notice however that this package was developed as part of our Master Thesis and because we have graduated and no longer work with this subject, it may take a while before issues are addressed. 

## References

Dua, D. & Graff, C. (2017). UCI machine learning repository.

Kovács, E. & Szántai, T. (2012). Hypergraphs as a mean of discovering the dependence structure of a discrete multivariate probability distribution. Ann Oper Res, 193, 71-90. <https://doi.org/10.1007/s10479-010-0814-y>

## See also

JOSS paper [![DOI](http://joss.theoj.org/papers/10.21105/joss.01480/status.svg)](https://doi.org/10.21105/joss.01480)
