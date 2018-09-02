#' Functions to solve problems of the form 
#' argmin sum <UX - Y> s.t. U^T G U = G,
#' where G is a self-adjoint matrix (G = G^T).
#' Typically, G = I_{p,q} and $U \in O(p,q)$, i.e.
#' the indefinite orthogonal group.
#' 
#' @author Joshua Agterberg
#' @references https://core.ac.uk/download/pdf/82794566.pdf and the dissertation referenced
#' therein.
library(Matrix)
library(alabama)


#' @param U in R^(n^2)
#' @param X the original X's
#' @param Y the original Y's
#' 
#' Function to get the value of f(U) and grad f(U).  X and Y are typically fixed
#' within the function.
f <- function(U,X,Y,A = NULL,b = NULL, gamma = NULL,grad = FALSE) {
  if ("Matrix" %in% class(U)) {
    stop('Error! U needs to be pre-vectorized (since it is symmetric)')
  }
  
  if (is.null(A)) {
    n <- dim(X)[2]
    A <- kronecker(Matrix::sparseMatrix(i= c(1:n),j=c(1:n),x=1),t(X%*%t(X)) )
  }
  
  if (is.null(b)) {
    b <- as.vector(Y %*% t(X))
  }
  
  if (is.null(gamma)) {
    gamma <- sum(diag(Y %*% t(Y)))
  }
  
  if (!grad) {
    return( sum( (A %*% U) * U) + - 2 * sum(b * U) + gamma)
  } else {
    #f_val <- sum( (A %*% U) * U) + - 2 * sum(b * U) + gamma
    grad_val <- 2*(A %*% U - b)
    return( as.vector(grad_val
    ))
  }
  
}

#' @param  U in R^(n^2)
#' Function to get the constraint vector and Matrix
constraint <- function(U,G,grad = FALSE) {
  U <- convertToMatrix(U)
  
  if (!grad) {
    
    
    return(as.vector(t(U) %*% G %*% U - G))
    
  } else {
    #cons <- as.vector(t(U) %*% G %*% U)
    grad_cons <- NULL
    k <- 1
    for (i in c(1:dim(U)[1])) {
      for (j in c(1:dim(U)[2])) {
        if (i == j) {
          J <- sparseMatrix(i = c(i),j = c(i),x=1,dims = dim(U))
        } else {
          J <- sparseMatrix(i = c(i,j),j = c(i,j),x=1,dims = dim(U))
        }
        
        grad_cons <- rbind(grad_cons,as.vector(G %*% U %*%  J))
        k <- k +1
      }
    }
    
    return(grad_cons)
  }
  
  
  
}

#main agorithm: method5.6
#iteration rule update: iterationRule 5.46 (eq)
#calculate dz: algorithm 5.8
# function to convert a symmetric matrix into an n(n+1)/2-length vector
# right now this basically does nothing special.
convertToVector <- function(B) {
  return(as.vector(B))
}

convertToMatrix<- function(b) {
  n <- length(b)
  return(Matrix(b,nrow = sqrt(n),ncol=sqrt(n)))
}


check <- function(n) {
  return(n*(n+1)/2)
}

solve_procrustes_problem <- function(X,Y,p,q) {
  solution <- list()
  X <- t(X)
  Y <- t(Y)
  
  #get the G matrix
  G <- Matrix::sparseMatrix(x=c(rep(1,p),rep(-1,q)),i = c(1:(p+q)),j = c(1:(p+q)))
  n <- dim(X)[1]
  A <- kronecker(Matrix::sparseMatrix(i= c(1:n),j=c(1:n),x=1),X%*%t(X))
  b <- as.vector(Y %*% t(X))
  gamma <- sum(diag(Y %*% t(Y)))
  f_new <- function(u) {
    return(f(u,X=X,Y=Y,A=A,b=b,gamma = gamma,grad=FALSE))
  }
  
  f_grad <- function(u) {
    return(f(u,X=X,Y=Y,A=A,b=b,gamma = gamma,grad=TRUE))
  }
  
  heq <- function(u) {
    return(constraint(u,G,grad = FALSE))
  }
  
  heq_jac <- function(u) {
    return(constraint(u,G,grad=TRUE))
  }
  
  
  I0 <- as.vector(sparseMatrix(i = c(1:dim(G)[1]),j=c(1:dim(G)[1]),x=1))
  
  
  result <- alabama::auglag(par = I0,fn = f_new,gr = f_grad,heq = heq, heq.jac = heq_jac)
  
  return(result) 
  
}

