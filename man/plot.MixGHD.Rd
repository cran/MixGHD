\name{plot}
\alias{plot.MixGHD}
\alias{plot,MixGHD-method}
\alias{plot,MixGHD,missing-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot objects of class MixGHD.
%%  ~~function to do ... ~~
}
\description{Plots the loglikelyhood function and for p<10 shows the splom of the data.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
	\S4method{plot}{MixGHD}(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A object of \code{\link{MixGHD-class}}};
  \item{y}{
Not used; for compatibility with generic plot.}
 %%     ~~Describe \code{x} here~~
}  
\details{ Plot the loglikhelyhood vale for each iteration of the EM algorithm. If p=2 it shows a contour plot. If 2<p<10 shows a splom of the data colored according to the cluster membership.
	}
\section{Methods}{
\describe{
\item{\code{signature(x = "MixGHD", y = "missing")}}{

S4 method for plotting objects of \code{\link{MixGHD-class}}.}
}
}
\author{Cristina Tortora.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{
 \code{\link{MixGHD-class}},\code{\link{MGHD}},\code{\link{MCGHD}},\code{\link{MSGHD}},\code{\link{cMSGHD}},\code{\link{MGHFA}}
}
\examples{
##loading banknote data
data(bankruptcy)


##model estimation
model=MCGHD(bankruptcy[,2:3],G=2,max.iter=30)

#result

plot(model)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Plot}

