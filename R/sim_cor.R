#' Simulation Cases
#'
#' \code{sim_cases} generates 7 cases with different covariance
#' correlation structures
#'
#' @param t is the case number
#' \enumerate{
#'   \item \eqn{t=1}: independent
#'   \item \eqn{t=2}: AR(1) with \eqn{\rho=0.5}
#'   \item \eqn{t=3}: AR(1) with \eqn{\rho=0.9}
#'   \item \eqn{t=4}: 2 compound symmetric blocks of size 3 with \eqn{\rho=0.5}
#'   \item \eqn{t=5}: 2 compound symmetric blocks of size 3 with \eqn{\rho=0.9}
#'   \item \eqn{t=6}: 3 compound symmetric blocks of size 2 with \eqn{\rho=0.5}
#'   \item \eqn{t=7}: 3 compound symmetric blocks of size 2 with \eqn{\rho=0.9}
#' }
#' @param n is the number of data points
#' @param b0 is the true coefficient vector
#' @param sd is the standard deviation
#' @export
sim_cases <- function(t,n=80,b0=c(1,1,1,0,0,0),sd=0.1){
  if(t==1){
    # independent
    data <- gen_exact(rho=0,p=6,k=6,n=n,bt=1,corstr="cs",sd=sd,b=b0)
  }
  if(t==2){
    # AR(1) rho=0.5
    data <- gen_exact(rho=0.5,p=6,k=6,n=n,bt=1,corstr="ar",sd=sd,b=b0)
  }
  if(t==3){
    # AR(1) rho=0.9
    data <- gen_exact(rho=0.9,p=6,k=6,n=n,bt=1,corstr="ar",sd=sd,b=b0)
  }
  if(t==4){
    # 2 compound symmetric blocks of size 3 rho=0.5
    data <- gen_exact(rho=0.5,p=6,k=3,n=n,bt=1,corstr="cs",sd=sd,b=b0)
  }
  if(t==5){
    # 2 compound symmetric blocks of size 3 rho=0.9
    data <- gen_exact(rho=0.9,p=6,k=3,n=n,bt=1,corstr="cs",sd=sd,b=b0)
  }
  if(t==6){
    # 3 compound symmetric blocks of size 2 rho=0.5
    data <- gen_exact(rho=0.5,p=6,k=2,n=n,bt=1,corstr="cs",sd=sd,b=b0)
  }
  if(t==7){
    # 3 compound symmetric blocks of size 2 rho=0.9
    data <- gen_exact(rho=0.9,p=6,k=2,n=n,bt=1,corstr="cs",sd=sd,b=b0)
  }
  return(data)
}

#' Generate simulations
#'
#' \code{gen_exact} is an internal function to \code{sim_cases} to
#' generate data with specific covariate correlations
#'
#' @param t is the case number
#' \enumerate{
#'   \item \eqn{t=1}: independent
#'   \item \eqn{t=2}: AR(1) with \eqn{\rho=0.5}
#'   \item \eqn{t=3}: AR(1) with \eqn{\rho=0.9}
#'   \item \eqn{t=4}: 2 compound symmetric blocks of size 3 with \eqn{\rho=0.5}
#'   \item \eqn{t=5}: 2 compound symmetric blocks of size 3 with \eqn{\rho=0.9}
#'   \item \eqn{t=6}: 3 compound symmetric blocks of size 2 with \eqn{\rho=0.5}
#'   \item \eqn{t=7}: 3 compound symmetric blocks of size 2 with \eqn{\rho=0.9}
#' }
#' @param rho is the correlation
#' @param p is the number of covariates
#' @param k is the number of covariates in a block
#' @param bt is the number of blocks
#' @param corstr is the correlation type ("ar" or "cs")
#' @param n is the number of data points
#' @param sd is the standard deviation
#' @param b is the true model coefficients
#' @export
gen_exact <- function(rho=0.7,p=2,k=2,bt=2,corstr="cs",n=100,sd=0.1,b=NA){
  # bt = 1: 1 pxp block
  # bt = 2: p/k blocks
  # bt = 3: 1 kxk block with the rest independent
  if(bt==1){
    if(corstr=="ar"){
      exponent <- abs(matrix(1:k - 1, nrow = k, ncol = k, byrow = TRUE) -
                        (1:k - 1))
      V0 <- rho^exponent
    }
    if(corstr=="cs"){
      V0 <- matrix(rho,ncol=k,nrow=k)
    }
    V <- kronecker(diag(p/k),V0)
    diag(V) <- 1
  }
  if(bt==2){
    V0 <- matrix(rho,ncol=k,nrow=k)
    V <- diag(p)
    V[1:k,1:k] <- V0
    diag(V) <- 1
  }

  ## Create independent covariates
  Xi <- kronecker(rep(1,10),hadamard(8)[,2:7])

  ## transform X
  X <- Xi%*%chol(V)

  Y <- X%*%b+rnorm(n,0,sd=sd)

  # appropriately scale data
  meanx <- colMeans(X)
  xc <- scale(X, meanx, FALSE)         # first subtracts mean
  normx <- apply(xc,2,norm,type="2")
  names(normx) <- NULL
  Xst <- scale(xc, FALSE, normx)
  meany <- mean(Y)
  Yst <- Y-meany
  #Yst <- (Y-meany)/(norm(Y,type="2"))
  return(list(X=X,Xst=Xst,Y=Y,Yst=Yst,beta_true=b))
}
