

dMSGHD <- function(data,p,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omegav=rep(1,p),lambdav=rep(0.5,p),gam=NULL,phi=NULL,log=FALSE) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  y=data
  ### y is the data set
  ##par is a list with all the parameter
  ## this is the density for 1 component so g is fixed
  # x is a n x p matrix
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  x     = y %*% (gam)  
  # mu    = par$mu; phi = par$phi;
  # alpha = par$alpha;
  alpha = alpha%*%gam
  mu =mu%*%gam
  d = p; chi = omegav; psi = omegav;
  lambda = lambdav;### getting the parameters from the input
  
  xmu = sweep(x,2,mu,"-")##x-mu
  
  pa = psi + alpha^2/phi   # p numeric
  mx = sweep(sweep(xmu^2,2,1/phi, FUN="*"), 2,chi, "+") # n x p matrix
  kx = sqrt(sweep(mx, 2, pa, "*")) # nxp matrix
  
  lx1 = sweep( sweep(log(mx),2,log(pa),"-"), 2, (lambda - 1/2)/2, "*")
  lx2 = t(apply(kx, 1, function(z,lam=NULL) { log(besselK( z, nu=lambda-1/2, expon.scaled =TRUE)) - z }, lam=lambda ))
  lx3 = sweep(xmu, 2, alpha/phi, FUN="*")
  
  lv  = matrix(0, nrow=d, ncol=3)###the density is divided in 3 parts
  lv[,1] = - 1/2*log( phi ) - 1/2*(log(2)+log(pi)) # d=1
  lv[,2] = - (log(besselK( sqrt(chi*psi), nu=lambda, expon.scaled =TRUE)) - sqrt(chi*psi) )
  lv[,3] = lambda/2*(log(psi)-log(chi) )
  
  if(ncol(y)==1){lx2=t(lx2)}
  val = apply(lx1 + lx2 + lx3, 1, sum) + sum(lv)
  
  if (!log) val = exp( val )
  
  return(val)
}


dCGHD <- function(data,p,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),lambda=1,omega=1,omegav=rep(1,p),lambdav=rep(1,p),wg=0.5,gam=NULL,phi=NULL) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
 
  

  wg=c(wg,1-wg)
  
  if (wg[2] != 0 & wg[1]!=0)  {
    val1 = dGHD(data,p,mu,alpha,sigma,omega,lambda)
    val2= dMSGHD(data,p,mu,alpha,sigma,omegav,lambdav,gam,phi)
    val =  wg[1]*(val1) + wg[2]*(val2)
  } else if(wg[1]==1) {val = dGHD(data,p,mu,alpha,sigma,omega,lambda)}
  else{val= dMSGHD(data,p,mu,alpha,sigma,omegav,lambdav,gam,phi)}
  
  return( val )
}




dGHD <- function(data, p, mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omega=1,lambda=0.5, log=FALSE) {
  # x  is a n * p matrix
  # generalized hyperbolic distribution
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  x=data
  
  # mu    = par$mu
  # sigma = par$sigma
  # alpha = par$alpha
  
  #check for input error
  if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
  d = length(mu)
  omega = exp(log(omega))

  
  
  #distances between points and mu
  invS = ginv(sigma)
  pa = omega + alpha %*% invS %*% alpha
  mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
  pa2=rep(pa,length(mx))
  kx = sqrt(mx*pa2)
  xmu = sweep(x,2,mu)
  
  lvx = matrix(0, nrow=nrow(x), 4)
  lvx[,1] = (lambda - d/2)*log(kx)
  lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
  lvx[,3] = xmu %*% invS %*% alpha
  
  
  lv = numeric(6)
  if(is.nan(log(det( sigma )))){sigma =diag(ncol(mu))}
  
  lv[1] = -1/2*log(det( sigma )) -d/2*(log(2)+log(pi))
  lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
  lv[3] = -lambda/2*( log(1) )
  lv[4] = lambda*log(omega) *0
  lv[5] = (d/2 - lambda)*log( pa )
  
  val = apply(lvx,1,sum) + sum(lv)
  if (!log) val = exp( val )
  
  return(val)
}

