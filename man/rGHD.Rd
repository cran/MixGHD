\name{rGHD}
\alias{rGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pseudo random number generation from a generalized hyperbolic distribution (GHD).
%%  ~~function to do ... ~~
}
\description{Generate n pseudo random numbers from a p dimensional generalized hyperbolic distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rGHD(n,p, mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),omega=1,lambda=0.5)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations.
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
}
\details{ The default values are: 0 for the mean and the skweness parameter alpha, diag(p) for sigma, 1 for omega, and 0.5 for lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A n times p matrix of numbers psudo randomly generated from a generilzed hyperbolic distribution}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\references{R.P. Browne, and  P.D. McNicholas (2015). A Mixture of  Generalized Hyperbolic Distributions. \emph{Canadian Journal of Statistics}, 43.2 176-198
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora, Aisha ElSherbiny, Ryan P. Browne, Brian C. Franczak, and Paul D. McNicholas.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{

data=rGHD(300,2,alpha=c(2,-2))

plot(data)
}
% R documentation directory.
\keyword{ Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

