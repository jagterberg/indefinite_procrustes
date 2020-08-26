library(ggplot2)
library(Rcpp)
Rcpp::sourceCpp("krein_inner_product.cpp")

krein_inner_prod <- function(x,y,p=length(x)-1) {
  return( kreinProduct(as.matrix(x),as.matrix(y),p)  )
}

krein_dist <- function(x,y,p=length(x)-1) {
  return( kreinProduct( as.matrix(x-y),as.matrix(x-y),p))
}

X0 <- c(.5,.3)
r <- .4
y1 <- seq(-1,1,.1)
y2 <- sort(rep(y1,length(y1)))
y <- rep(y1,length(y2))
yn <- cbind(y,y2)

all_probs <- kreinProduct(yn,yn,1)


ball_vals <- apply(yn,1,function(y) {
    return(krein_dist(X0,y))
})

X0_new <- matrix(0,nrow(yn),ncol(yn))
X0_new[,1] <- X0[1]
X0_new[,2] <- X0[2]

ball_vals <- krein_dist(X0,yn)
ball_probs <- krein_inner_prod(yn,yn)
inds <- (ball_probs <= 1 & ball_probs >= 0)
ball_vals <- ball_vals[inds]
yn <- yn[inds,]
shade <- ifelse(ball_vals < r^2 & ball_vals > -(r^2),1,0)
dat <- data.frame(y1 = yn[,1],y2 = yn[,2],shade=factor(shade))
g <- ggplot(data =dat ) + geom_point(aes(x=y1,y=y2,col=shade))#,alpha=shade2))
g


  
           
