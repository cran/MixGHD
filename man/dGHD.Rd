\name{dGHD}
\alias{dGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density of a generalized hyperbolic distribution (GHD).
%%  ~~function to do ... ~~
}
\description{Compute the density of a p dimensional generalized hyperbolic distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dGHD(data,p, mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omega=1,lambda=0.5, log=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ n x p data set
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
 \item{omega}{ (optional) the unidimensional concentration parameter omega
%%     ~~Describe \code{q} here~~
}
  \item{lambda}{ (optional) the unidimensional index parameter lambda
%%     ~~Describe \code{q} here~~
}
  \item{log}{ (optional) if TRUE returns the log of the density
%%     ~~Describe \code{q} here~~
}
}
\details{ The default values are: 0 for the mean and the skweness parameter alpha, diag(p) for sigma, 1 for omega, and 0.5 for lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A n dimensional vector with the density from   a generilzed hyperbolic distribution}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{ R.P. Browne, and  P.D. McNicholas (2015). A Mixture of  Generalized Hyperbolic Distributions. \emph{Canadian Journal of Statistics}
%% ~put references to the literature/web site here ~
}
\author{Ryan P. Browne, Cristina Tortora
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{



x = seq(-3,3,length.out=50)
y = seq(-3,3,length.out=50)
xyS1 = matrix(0,nrow=length(x),ncol=length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
      xy <- matrix(cbind(x[i],y[j]),1,2)	
      xyS1[i,j] =  dGHD(xy,2) 
      
    }
  }
contour(x=x,y=y,z=xyS1, levels=c(.005,.01,.025,.05, .1,.25), main="MGHD",ylim=c(-3,3), xlim=c(-3,3))




}
% R documentation directory.
\keyword{ Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

