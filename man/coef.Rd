\name{coef}
\alias{coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coefficients
%%  ~~function to do ... ~~
}
\description{Coefficents of the estimated model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coef(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{input}{ An S4 object of class \link{MixGHD}.%%     ~~Describe \code{x} here~~
}
}
\value{The coefficents of the estimated model
	}
\author{Cristina Tortora
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{
##loading bankruptcy data
data(bankruptcy)

##model estimation
res=MCGHD(data=bankruptcy[,2:3],G=2,method="kmedoids",max.iter=30)
#rcoefficients of the model
coef(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

