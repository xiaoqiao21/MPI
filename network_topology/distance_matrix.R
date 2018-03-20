library(igraph)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('network topology in super computer.cpp')
a1 <- a31(64)
g31 <- graph_from_adjacency_matrix(a1)
disa31 <- distances(g31)
eigen31 <- eigen(a1)$va
a2 <- a32(8,64)
g32 <- graph_from_adjacency_matrix(a2)
disa32 <- distances(g32)
eigen32 <- eigen(a2)$va
a3 <- a33(12)
g33 <- graph_from_adjacency_matrix(a3)
disa33 <- distances(g33)
eigen33 <- eigen(a3)$va