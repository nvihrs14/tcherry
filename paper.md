---
title: 'tcherry: Learning the structure of tcherry trees'
authors:
- affiliation: 1
  name: Katrine Kirkeby
  orcid: 0000-0001-6065-1583
- affiliation: 1
  name: Maria Knudsen
  orcid: 0000-0002-0067-6671
- affiliation: 1
  name: Ninna Vihrs
  orcid: 0000-0002-7507-8473
date: "26 April 2019"
bibliography: paper.bib
tags:
- t-cherry trees
- Undirected triangulated graphical models
- Dependencies between variables
affiliations:
- index: 1
  name: Department of Mathematical Sciences, Aalborg University
---

# Summary

The R [@R] package `tcherry` contains a variety of functions for learning the structure of a k'th order t-cherry tree from given categorical data, see for instance @EKTShyp for an explanation of this concept. This is a graphical model extending what is known as a Chow-Liu tree [@ChowLiu]. Chow-Liu trees have for instance been used to estimate population frequencies of Y-STR haplotypes in @Y-STR and t-cherry trees have been used to model relationships in social networks in @Proulx. The functions attempt to find a t-cherry tree structure of maximal likelihood. To do this exact, it is necesarry to investigate all possible t-cherry tree structures of the given order. This is in most cases to time-consuming and therefore most of the functions use greedy search algorithms. Some implementations are inspired by algorithms in @EKTS, @EKTSdisc and @Proulx, but the package also contains some new algorithms and extensions. The package is only for structure learning and only categorical data is supported. 

The package can be used as a tool to analyse problems exploring dependencies between any kind of categorical variables. The fitted t-cherry structure can be used to make statements about conditional dependencies and independencies. The structure can also be used for pattern recognition and independence statements can be used for variabel selection for a prediction problem [@patternrec]. If the structure is used in combination with packages such as `gRain` [@gRain], it may also be used to estimate probability distributions of the variables or for prediction. This makes it possible to use the structure as an expert system.

The t-cherry tree structure can be used in a variety of scientific fields such as biostatistics and artificial intelligence. The audience of the package is anyone who wants to model dependencies between categorical variables, approximate their probability distribution or solve klassification problems with categorical variables.

The following figure shows an example of a fourth order t-cherry tree learned from the car evaluation data set from UCI Machine Learning Repository [@UCI].

![](inst/images/tch.png)

The R package `tcherry` is available on [GitHub](https://github.com/nvihrs14/tcherry).

# References
