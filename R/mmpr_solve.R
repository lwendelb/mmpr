#' Ridge Regression
#'
#' \code{z_fun} is an internal function to \code{coord_desc}. It returns
#' the normalizing constant for a covariate
#' (returns 1 if X is normalized to have mean 0
#' and 2-Norm 1)
#'
#'
#' @param x design matrix
#' @param ind scalar covariate index
#' @export
z_fun <- function(x,ind){
  sum(x[,ind]^2)
}

#' \code{rho_fun} is an internal function to \code{coord_desc}. It returns
#' the partial sum of squared errors
#'
#' @param x design matrix
#' @param y response vector
#' @param b vector of coefficients for a single model
#' @param ind scalar covariate index
#' @export
rho_fun <- function(x,y,b,ind){
  val <- t(x[,ind])%*%(y-x[,-ind]%*%b[-ind])
  return(val)
}

#' \code{gamma_fun} is an internal function to \code{coord_desc}. It returns
#' the value to softthreshold by
#'
#' @param bp \eqn{(p \times (M-1)} matrix of models
#' @param lambda scalar sparsity penalty weight
#' @param w scalar similarity penalty weight
#' @param ind scalar covariate index
#' @param M number of models
#' @param p number of covariates
#' @param c scalar sparsity setting, \eqn{c=0,1}
#' @param d scalar similarity setting, \eqn{d=0,1}
#' @export
gamma_fun <- function(bp,lambda,w,ind,M,p,c,d){
  bp <- matrix(bp,nrow=(M-1),ncol=p)
  val <- lambda*(2-c)+(2-d)*w*sum(abs(bp[,ind])^d)
  return(val)
}

#' \code{theta_fun} is an internal function to \code{coord_desc}. It returns
#' the value to shrink by
#'
#' @param bp \eqn{(p \times (M-1)} matrix of models
#' @param lambda scalar sparsity penalty weight
#' @param w scalar similarity penalty weight
#' @param ind scalar covariate index
#' @param z normalizing constant
#' @param M number of models
#' @param p number of covariates
#' @param c scalar sparsity setting, \eqn{c=0,1}
#' @param d scalar similarity setting, \eqn{d=0,1}
#' @export
theta_fun <- function(bp,lambda,w,ind,z,M,p,c,d){
  bp <- matrix(bp,nrow=(M-1),ncol=p)
  val <- z+lambda*(c-1)+(d-1)*w*sum(abs(bp[,ind])^d)
  return(val)
}

#' \code{Soft} is an internal function to \code{coord_desc}. It returns
#' the shrunken soft threshold of rho
#'
#' @param g scalar threshold
#' @param t scalar shrinkage
#' @param rho value to perform soft-thresholding on
#' @export
Soft <- function(g,t,rho){
  (1/t)*sign(rho)*max(0,(abs(rho)-g/2))
}

#' \code{multi_model} is a general function to perform mmpr
#'
#' @param X \eqn{(n \times p} covariate matrix
#' @param Y \eqn{(n \times 1)} response
#' @param w scalar similarity penalty weight
#' @param lambda scalar sparsity penalty weight
#' @param init initial model
#' @param M number of models
#' @params steps maximum number of iterations to perform
#' @param zero_thresh number of steps after which a coefficient
#' equal to zero stops changing value
#' @params tol tolerance
#' @param c scalar sparsity setting, \eqn{c=0,1}
#' @param d scalar similarity setting, \eqn{d=0,1}
#' @export
multi_model <- function(X,Y,w=0,lambda=0,init=0,M=2,method="coord_desc",steps=5000, zero_thresh = 50,tol = 1e-4,c=NA,d=NA,iter=NA){
  if(method=="coord_desc"){
    return(coord_desc(X,Y,w=0,lambda=0,init=0,M=2,steps=5000, zero_thresh = 50,tol = 1e-4,c=NA,d=NA,iter=NA))
  }
  else{
    stop("specify method")
  }
}

#' \code{coord_desc} performs mmpr via cyclic coordinate descent
#'
#' @inheritParams multi_model
#' @export
coord_desc <- function(X,Y,w=0,lambda=0,init=0,M=2,steps=5000, zero_thresh = 50,tol = 1e-4,c=NA,d=NA,iter=NA){
  if(is.na(c)|is.na(d)){
    c=1
    d=1
    warning("c and d unspecified in w_opt, defaulting to c=1, d=1")
  }
  if(w<0){
    stop("tuning parameter w is negative")
  }
  if(lambda<0){
    stop("tuning parameter w is negative")
  }
  p <- ncol(X)
  b_ca <- array(0,dim=c(p,M,steps))
  b_ca[,,1] <- init
  for(step in 2:steps){
    b_ca[,,step] <- b_ca[,,(step-1)] # update current candidates with last round
    for(i in 1:p){ # cycle through parameters in each model
      for(j in 1:M){ # cycle through models

        ## Calculate necessary quantites
        beta <- b_ca[,,step]
        if(b_ca[i,j,step]==0 & (step > zero_thresh)){
          b_ca[i,j,step] <- 0
        }
        else{
          zj <- z_fun(x=X,ind=j)
          rhoj <- rho_fun(x=X,y=Y,b=beta[,i],ind=j)
          gammaj <- gamma_fun(bp=beta[,-i],lambda=lambda,w=w,ind=j,M=M,p=p,c=c,d=d)
          thetaj <- theta_fun(bp=beta[,-i],lambda=lambda,w=w,ind=j,z=zj,M=M,p=p,c=c,d=d)
          bhat <- Soft(g=gammaj,t=thetaj,rho=rhoj)
          b_ca[i,j,step] <- bhat
        }
      }
    }
    #print(max(abs(b_ca[,,step]-b_ca[,,(step-1)])))
    if(all(abs(b_ca[,,step]-b_ca[,,(step-1)])<tol)){
      return(b_ca[,,step])
    }
  }
  warning(paste("max iters reached for iter",iter,", variable changing by ",max(abs(b_ca[,,step]-b_ca[,,(step-1)]))))
  return(b_ca[,,steps])
}
