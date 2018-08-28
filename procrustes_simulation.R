source("C:/Users/joshu/Documents/Research/procrustes/indefinite_procrustes.R")

library(irlba)

simulateAdjacency <- function(B = B_mat,n = 100) {
  A <- matrix(0,n,n)
  for (i in c(2:n)) {
    for (j in c(i:(n-1))) {
      if(i < n/2) {
        if (j < n/2) {
          A[i,j] <- rbinom(1,1,B[1,1])
        } else {
          A[i,j] <- rbinom(1,1,prob=B[1,2])
        }
      } else {
        if( j > n/2) {
          A[i,j] <- rbinom(1,1,prob=B[2,2])
        } else {
          A[i,j] <- rbinom(1,1,prob=B[1,2])
        }
      }
      A[j,i] <- A[i,j]
    }
  }
  
  return(A)
}

procrustes_plot <- function(Xhat,X_true,X_new) {
  n <- dim(Xhat)[1]
  dat <- data.frame(X1 = c(X_true[,1],X_new[,1],Xhat[,1])
                    ,X2 =c(X_true[,2],X_new[,2],Xhat[,2])
                    ,Type = c(rep('True',n),rep('Rotated',n),rep('Original',n))
                    ,Indicator = c(rep(1,n),rep(0,2*n)))
  
  g <- ggplot(data=dat) + geom_point(aes(x = X1,y=X2,col = Type,size = factor(Indicator)
                                         ,alpha = factor(Indicator))) + 
    scale_color_manual(values = c( 'Rotated' = 'red','Original' ='orange',
                                   'True' = 'black'),
                       labels= c('Rotated' ='Xhat After Transform'
                                 ,'Original' = ' Xhat Original',
                                 'True' = 'True Values')
    ) + guides(colour = guide_legend('')) +
    scale_alpha_manual(values = c(.3,1),guide=FALSE) +
    scale_size_manual(values = c('1' = 4,'0'= 2),guide = FALSE) +
    ggtitle('Spectral Embedding of Indefinite Matrix after
            Applying Indefinite Procrustes Analysis')
  
  print(g)
}

get_true_X <- function(n,B = B_mat) {
  X <- matrix(0,n,n)
  X[c(1:(n/2)),c(1:(n/2))] <- B[1,1]
  X[c(1:(n/2)),c(((n/2) + 1):n)] <- B[1,2]
  X[c(((n/2)+1):n),c(1:(n/2))] <- B[2,1]
  X[c(((n/2)+1):n),c(((n/2)+1):n)] <- B[2,2]
  
  X_true <- irlba(X,2)
  X_true <- X_true$u[,c(1:2)] %*% diag(X_true$d[c(1:2)])^(.5)
  return(X_true)
}
#simulate one matrix:
p <- .2; q <- .7; r<- .2
B_mat <- matrix(c(p,q,q,r),2,2)
p <- 1
q <- 1
n <- 100
A <- simulateAdjacency(n=n)
X_true <- get_true_X(n)



A_svd <- irlba(A,2)
Xhat <- A_svd$u[,c(1:2)] %*% diag(A_svd$d[c(1:2)])^(.5)
X <- Xhat

Y <- X_true
solution <- solve_procrustes_problem(X,Y,p=1,q=1)

X1 <- X %*% convertToMatrix(solution$par)
X2 <- Y


library(ggplot2)
procrustes_plot(Xhat = Xhat,X_true = X2,X_new = X1)

n <- 100
A1 <- simulateAdjacency(n=n)
A2 <- simulateAdjacency(n=n)
A1_svd <- irlba(A1,2)
A2_svd <- irlba(A2,2)
X1 <- A1_svd$u[,c(1:2)] %*% diag(A1_svd$d[c(1:2)])^(.5)
X2 <- A2_svd$u[,c(1:2)] %*% diag(A2_svd$d[c(1:2)])^(.5)

X_true <- get_true_X(n)

solution1 <- solve_procrustes_problem(X1,X2,p=1,q=1)
solution2 <- solve_procrustes_problem(X1,X2,p=2,q=0)

procrustes_plot(Xhat = (X1 %*% convertToMatrix(solution1$par)),
                X_true= X_true, X_new = X2)
procrustes_plot(Xhat = (X1 %*% convertToMatrix(solution2$par)),
                                X_true= X_true, X_new = X2)




