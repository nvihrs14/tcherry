cond_independence_test <- function(var1, var2, cond = c(), data,
                                   smooth = 0){
  if (length(cond) == 0){
    ncond <- 1

    tab_12 <- table(data[, c(var1, var2)])
    p_tab_12 <- (tab_12 + smooth) / sum(tab_12 + smooth)

    tab_1 <- table(data[, var1])
    names(dimnames(tab_1)) <- var1
    p_tab_1 <- (tab_1 + smooth) / sum(tab_1 + smooth)

    tab_2 <- table(data[, var2])
    names(dimnames(tab_2)) <- var2
    p_tab_2 <- (tab_2 + smooth) / sum(tab_2 + smooth)

    chi_sq <- 2 *
      sum(gRbase::ar_prod(tab_12,
                          log(gRbase::ar_div(p_tab_12,
                                             gRbase::ar_prod(p_tab_1,
                                                             p_tab_2)))))
  }else{
    tab_cond <- table(data[, cond], dnn = cond)
    p_tab_cond <- (tab_cond + smooth) / sum(tab_cond + smooth)

    ncond <- length(tab_cond)

    tab_12cond <- table(data[, c(var1, var2, cond)])
    p_tab_12cond <- (tab_12cond + smooth) / sum(tab_12cond + smooth)
    p_tab_12_given_cond <- gRbase::ar_div(p_tab_12cond, p_tab_cond)

    tab_1cond <- table(data[, c(var1, cond)])
    p_tab_1cond <- (tab_1cond + smooth) / sum(tab_1cond + smooth)
    p_tab_1_given_cond <- gRbase::ar_div(p_tab_1cond, p_tab_cond)

    tab_2cond <- table(data[, c(var2, cond)])
    p_tab_2cond <- (tab_2cond + smooth) / sum(tab_2cond + smooth)
    p_tab_2_given_cond <- gRbase::ar_div(p_tab_2cond, p_tab_cond)

    chi_sq <- 2 *
      sum(gRbase::ar_prod(tab_12cond,
                          log(gRbase::ar_div(
                            p_tab_12_given_cond,
                            gRbase::ar_prod(p_tab_1_given_cond,
                                            p_tab_2_given_cond)))))
  }

  n1 <- length(unique(data[, var1]))
  n2 <- length(unique(data[, var2]))

  df <- (n1 - 1) * (n2 - 1) * ncond

  p <- 1 - pchisq(chi_sq, df)

  return(list("chi_sq_statistic" = chi_sq,
              "df" = df,
              "p_value" = p))

}
