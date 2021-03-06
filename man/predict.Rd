\name{predict}
\alias{predict.MixGHD}
\alias{predict,MixGHD-method}
\alias{predict,MixGHD,missing-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Membership prediction for objects of class MixGHD
%%  ~~function to do ... ~~
}
\description{Cluster membership 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
	\S4method{predict}{MixGHD}(object)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{ An S4 object of class \link{MixGHD}.%%     ~~Describe \code{x} here~~
}
}
\value{The cluster membership 
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
predict(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

