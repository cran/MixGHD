\name{dCGHD}
\alias{dCGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of a coalesced generalized hyperbolic distribution (MSGHD).
%%  ~~function to do ... ~~
}
\description{Compute the density of a p dimensional  coalesced generalized hyperbolic distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dCGHD(data,p,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),lambda=1,omega=1,
  omegav=rep(1,p),lambdav=rep(1,p),wg=0.5,gam=NULL,phi=NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{n x p data set
%%     ~~Describe \code{x} here~~
}
  \item{p}{ number of variables.
%%     ~~Describe \code{x} here~~
}
  \item{mu}{ (optional) the p dimensional mean
%%     ~~Describe \code{q} here~~
}
 \item{alpha}{ (optional) the p dimensional skewness parameter alpha
%%     ~~Describe \code{q} here~~
}
 \item{sigma}{ (optional) the p x p dimensional scale matrix
%%     ~~Describe \code{q} here~~
}
  \item{lambda}{ (optional) the 1 dimensional index parameter lambda
%%     ~~Describe \code{q} here~~
}
\item{omega}{ (optional) the 1 dimensional concentration parameter omega
%%     ~~Describe \code{q} here~~
}
 \item{omegav}{ (optional) the p dimensional concentration parameter omega
%%     ~~Describe \code{q} here~~
}
  \item{lambdav}{ (optional) the p dimensional index parameter lambda
%%     ~~Describe \code{q} here~~
}
 \item{wg}{ (optional) weight
%%     ~~Describe \code{q} here~~
}
  \item{gam}{ (optional) the pxp gamma matrix
%%     ~~Describe \code{q} here~~
}
  \item{phi}{ (optional) the p dimensional vector phi
%%     ~~Describe \code{q} here~~
}
}
\details{ The default values are: 0 for the mean and the skweness parameter alpha, diag(p) for sigma, 1 for omega, and 0.5 for lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A n dimensional vector with the density from a  coalesced generilzed hyperbolic distribution}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{C. Tortora, B.C. Franczak, R.P. Browne, A. ElSherbiny, and  P.D. McNicholas (2014). Mixture of Multiple Scaled Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1403.2332v2
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Brian C. Franczak 
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{



x = seq(-3,3,length.out=30)
y = seq(-3,3,length.out=30)
xyS1 = matrix(0,nrow=length(x),ncol=length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
      xy <- matrix(cbind(x[i],y[j]),1,2)	
      xyS1[i,j] =  dCGHD(xy,2) 
      
    }
  }
contour(x=x,y=y,z=xyS1, levels=c(.005,.01,.025,.05, .1,.25), main="CGHD",ylim=c(-3,3), xlim=c(-3,3))


}

% R documentation directory.
\keyword{ Coalesced Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

