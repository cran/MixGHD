rGHD <- function(n,p, mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omega=1,lambda=0.5) {

  d = p
  sigma = sigma
  alpha = alpha
  mu = mu
  chi = omega
  psi =omega
  lambda = lambda
  
  #print(c(chi,psi,lambda))
  
  w = rgig(n = n, lambda = lambda, chi = chi, psi = psi)
  x = rmvnorm(n, sigma=(t(sigma)%*%sigma))
  x = sweep(x, 1, sqrt(w), FUN="*")
  
  A = matrix( alpha, nrow=n, ncol=length(alpha), byrow=TRUE)
  A = sweep(A, 1, w, FUN="*")
  
  
  y = A + x
  y = sweep(y, 2, mu, FUN="+") 
  val = y#list(x=y,w=w)
  return(val)
}

##MSGH
rMSGHD <- function(n,p,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omegav=rep(1,p),lambdav=rep(0.5,p)) {
  
  p = p
  x = matrix(0, nrow=n, ncol=p)
  y = matrix(0, nrow=n, ncol=p)
  sigma = (sigma)
  alpha = alpha
  mu = mu
  chi = omegav
  psi = omegav
  lambda = lambdav
  phi=diag(p)
  diag(phi)=(eigen( sigma)$values)
  gam=eigen( sigma)$vectors
  #print(c(chi,psi,lambda))
  for (i in 1:n) {
    ww=matrix(0,p,1)
    for(j in 1:p){
      w = rgig(n = 1, lambda = lambda[j], chi = chi[j], psi = psi[j])
      ww[j,1]=w}
    wd=diag(p)
    diag(wd)=ww
    x[i,] = t(wd%*%(alpha))+t(sqrt(wd)%*%t((rnorm(p))%*%(gam%*%sigma%*%phi%*%t(gam))))
  }
  
  
  xn=x+t(replicate(n,mu))
  y = xn %*% t(gam)
  val = xn#list(y=y, x=xn,w=w)
  return(val)
}

rCGHD <- function(n,p,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omega=1,lambda=0.5,omegav=rep(1,p),lambdav=rep(0.5,p),wg=0.5) {
  gh=rGHD(n,p, mu=mu,alpha=alpha,sigma=sigma,omega=omega,lambda=lambda) 
  ms=rMSGHD(n,p, mu=mu,alpha=alpha,sigma=sigma,omegav=omegav,lambdav=lambdav) 
  cgh=wg*gh+(1-wg)*ms
  x=cgh
  return(x)
}