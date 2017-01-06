\name{contourpl}
\alias{contourpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Contour plot
%%  ~~function to do ... ~~
}
\description{Contour plot for a given set of parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
contourpl(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{input}{ An S4 object of class \link{MixGHD}.%%     ~~Describe \code{x} here~~
}
}
\value{The contour plot
	}
\author{Cristina Tortora
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{
##loading bankruptcy data
data(bankruptcy)

##model estimation
res=MCGHD(data=bankruptcy[,2:3],G=2,method="kmedoids",max.iter=30)
#result
contourpl(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

