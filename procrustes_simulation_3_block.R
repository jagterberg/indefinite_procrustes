#' ---
#' title: "Indefinite Procrustes"
#' author: "Joshua Agterberg"
#' date: ""
#' 
#'output: html_file
#' ---
#' 
#' 
#' 
#' 
#' First, load dependencies:

source("C:/Users/joshu/Documents/Research/procrustes/indefinite_procrustes.R")
library(Matrix)
library(irlba)
library(ggplot2)


simulateAdjacency <- function(B = B_mat,n = 99,b= c(rep(1,n/3),rep(2,n/3),rep(3,n/3))) {
  A <- matrix(0,n,n)
  for (i in c(1:n)) {
    j <- n
    while (j > i) {
      A[i,j] <- rbinom(n=1,size=1,p=B[b[i],b[j]])
      A[j,i] <- A[i,j]
      j <- j -1
    }
  }
  
  return(A)
}


#' Function to plot X, Y, and some truth value.  Assumed to be the first
#' two columns of X and Y, although others can be provided through which.inds.
#' Also, names are assumed to be the in the form( Rotated, Original,True), but
#' the names can be changed to different values through the names parameter.
#' Finally, the truth.ind parameter specifies a length 3n vector of an indicator
#' for the truth -- this will be the size of the truth value in the plot.
procrustes_plot <- function(True,X,Y,which.inds = c(1,2),names = NULL,truth.ind = NULL
                            , title = c('Spectral Embedding of Indefinite Matrix after
                                        Applying Indefinite Procrustes Analysis')
                            ,transparency = .2,true.size = 4) {
  n <- dim(X)[1]
  if (is.null(names)) {
    names <- factor(c(rep('Rotated',n),rep('Original',n),rep('True',n)),
                    levels=c('Rotated','Original','True'))
  }
  
  if (is.null(truth.ind)) {
    truth.ind <- factor(c(rep(0,2*n),rep(1,n)),levels = c("0","1"))
  }
  
  dat <- data.frame(X1 = c(X[,which.inds[1]],Y[,which.inds[1]],True[,which.inds[1]])
                    ,X2 =c(X[,which.inds[2]],Y[,which.inds[2]],True[,which.inds[2]])
                    ,Type = names
                    ,Indicator =truth.ind )
  names.levels <- levels(names)
  truth.ind.levels <- levels(truth.ind)
  truths.size<- c(1.5,true.size)
  truths.alpha <- c(transparency,1)
  names(truths.alpha) <- truth.ind.levels
  names(truths.size) <- truth.ind.levels
  names.cols <- c('red','orange','black')
  names(names.cols) <- names.levels
  names.labels <- c('Xhat after Rotation','Xhat Original','True Values')
  names(names.labels) <- names.levels
  g <- ggplot(data=dat) + geom_point(aes(x = X1,y=X2,col = Type,size = as.factor(Indicator)
                                         ,alpha = factor(Indicator))) + 
    scale_color_manual(values = names.cols,
                       labels= names.labels
    ) + guides(colour = guide_legend('')) +
    scale_alpha_manual(values = truths.alpha,guide=FALSE) +
    scale_size_manual(values = truths.size
                      ,guide = FALSE) +
    ggtitle(title)
  
  print(g)
}


#' Function to create the true latent position vectors from an SBM 
#' with block probabilities B and latent assignment vector b
get_true_X <- function(n,B,b= c(rep(1,n/3),rep(2,n/3),rep(3,n/3))) {
  X <- matrix(0,n,n)
  for (i in c(1:n)) {
    j <- n
    while (j >= i) {
      X[i,j] <- B[b[i],b[j]]
      X[j,i] <- X[i,j]
      j <- j - 1
    }
  }
  d <- length(unique(b))
  X_true <- irlba(X,d)
  X_true <- X_true$u[,c(1:d)] %*% diag(X_true$d[c(1:d)])^(.5)
  return(X_true)
}


#' Function to plot procrustes for two block model:
procrustes.simulation2 <- function(seed,n,p = .2,q=.7,r=.2) {
  set.seed(seed)
  B_mat <- matrix(c(p,q,q,r),2,2)
  b <- c(rep(1,n/2),rep(2,n/2))
  A <- simulateAdjacency(n=n,B= B_mat,b=b)
  Xtrue <- get_true_X(n,B = B_mat,b=b)
  A_svd <- irlba(A,2)
  Xhat <- A_svd$u[,c(1:2)] %*% diag(A_svd$d[c(1:2)])^(.5)
  solution <- solve_procrustes_problem(X=Xhat,Y=Xtrue,p=1,q=1)
  
  Xhat_rotated <- Xhat %*% convertToMatrix(solution$par)
  
  procrustes_plot(True = Xtrue,Y = Xhat,X = Xhat_rotated,transparency = .3,true.size= 3)
}

procrustes.simulation2(24,1000)


procrustes.simulation3 <- function(seed,n
                                   ,B_mat =  matrix(c(.6,.9,.9,.9,.6,.9,.9,.9,.3),3,3)
                                   ,dims = c(1,2)) {
  set.seed(seed)
  b <- c(rep(1,n/3),rep(2,n/3),rep(3,n/3))
  A <- simulateAdjacency(n=n,B= B_mat,b=b)
  Xtrue <- get_true_X(n,B = B_mat,b=b)
  A_svd <- irlba(A,3)
  Xhat <- A_svd$u[,c(1:3)] %*% diag(A_svd$d[c(1:3)])^(.5)
  solution <- solve_procrustes_problem(X=Xhat,Y=Xtrue,p=1,q=2)
  
  Xhat_rotated <- Xhat %*% convertToMatrix(solution$par)
  title <- paste0("Procrustes Analysis For D=3 \n(on x = dimension ",dims[1],", y = dimension ",dims[2],")")
  procrustes_plot(True = Xtrue[,dims],Y = Xhat[,dims],X = Xhat_rotated[,dims]
                  ,transparency = .3,true.size= 3
                  ,title =title)
  
}
seed = 100
procrustes.simulation3(seed=seed,n=666,dims = c(1,2))
procrustes.simulation3(seed=seed,n=666,dims = c(1,3))
procrustes.simulation3(seed=seed,n=666,dims = c(2,3))


seed =300

procrustes.simulation3(seed=seed,n=666,dims = c(1,2))
procrustes.simulation3(seed=seed,n=666,dims = c(1,3))
procrustes.simulation3(seed=seed,n=666,dims = c(2,3))


library(rgl)
procrustes_plot3D <- function(True,X,Y,names = NULL,truth.ind = NULL
                            , title = c('Spectral Embedding of Indefinite Matrix after
                                        Applying Indefinite Procrustes Analysis')
                            ,transparency = .2,true.size = 4) {
  n <- dim(X)[1]
  if (is.null(names)) {
    names <- factor(c(rep('Rotated',n),rep('Original',n),rep('True',n)),
                    levels=c('Rotated','Original','True'))
  }
  
  if (is.null(truth.ind)) {
    truth.ind <- factor(c(rep(0,2*n),rep(1,n)),levels = c("0","1"))
  }
  
  dat <- data.frame(X1 = c(X[,1],Y[,1],True[,1])
                    ,X2 =c(X[,2],Y[,2],True[,2])
                    ,X3 = c(X[,3],Y[,3],True[,3])
                    ,Type = names
                    ,Indicator =truth.ind )
 
  dat$Type2 <- ifelse(dat$Type == "Rotated","red",
                     ifelse(dat$Type == "True","black","lightblue"))
  title <- "Indefinite Procrustes Analysis From GRDPG 3-block Model"
  title2 <- "lightblue: Non-Rotated, Red: Rotated, Black: True"
  with(dat,plot3d(X1,X2,X3,col = dat$Type2,size = 0,main=title,sub=title2,box=FALSE,axes=TRUE))
 
  irisList <- split(dat,dat$Indicator)
  
  # Setup the plot

  # Use a separate call to points3d() to plot points of each size
  for(i in seq_along(irisList)) {
    if (irisList[[i]]$Indicator[1] == "1") {
      j <- 10
    
    } else {
      j <- 2
    }
    with(irisList[[i]], points3d(X1,X2,X3, col=Type2, size=j))
  }
  

  
}

procrustes.simulation3D <- function(seed,n
                                    ,B_mat =  matrix(c(.6,.9,.9,.9,.6,.9,.9,.9,.3),3,3)
                                    ) {
  set.seed(seed)
  b <- c(rep(1,n/3),rep(2,n/3),rep(3,n/3))
  A <- simulateAdjacency(n=n,B= B_mat,b=b)
  Xtrue <- get_true_X(n,B = B_mat,b=b)
  A_svd <- irlba(A,3)
  Xhat <- A_svd$u[,c(1:3)] %*% diag(A_svd$d[c(1:3)])^(.5)
  solution <- solve_procrustes_problem(X=Xhat,Y=Xtrue,p=1,q=2)

  Xhat_rotated <- Xhat %*% convertToMatrix(solution$par)
  
  #title <- paste0("Procrustes Analysis For D=3 \n(on x = dimension ",dims[1],", y = dimension ",dims[2],")")
  
  procrustes_plot3D(True = Xtrue,Y = Xhat,X = Xhat_rotated
                    ,transparency = .3,true.size= 3
                    ,title =title)
}

procrustes.simulation3D(100,666)

