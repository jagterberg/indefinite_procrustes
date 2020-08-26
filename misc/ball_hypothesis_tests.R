library(Rcpp)
library(Matrix)
library(irlba)
Rcpp::sourceCpp("C:/Users/joshu/Documents/Research/generate_adjacency_matrix.cpp")
source("C:/Users/joshu/Documents/Research/procrustes/indefinite_procrustes.R")

n <- 1000
B <- matrix(c(.4,.1,.1,.5),2,2)
b <- c(rep(1,n/2),rep(2,n/2))
P <- getProbabilityMatrix(b,n,B)
true_x <- irlba(P,2)
true_x <- true_x$u[,1:2] %*% diag((true_x$d)^.5)[,1:2]
d <- 2

reps <- 100
set.seed(1234)
testStats <- rep(0,reps)
for (k in c(1:reps)) {
  A <- generateAdjacencyMatrix(P)
  X <- irlba(A,2)
  Xhat <- X$u[,1:2] %*% diag((X$d)^.5)[,1:2]
  M <- t(true_x) %*% Xhat
  gamma2 <- X$d[2]
  M <- svd(M)
  O <- M$u %*% t(M$v)
  testStat <- sum(
    (Xhat %*% O - true_x)^2
  )
  
  testStats[k] <- testStat
}

mean(testStats)
