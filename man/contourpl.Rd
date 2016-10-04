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
contourpl(data,output,method="MGHD",scale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{ A n x p matrix or data frame such that rows correspond to observations and columns correspond to variables, the function only uses the first 2 variables.
%%     ~~Describe \code{x} here~~
}
  \item{output}{  output of one of the functions MGHD, MCGHD, MSGHD, cMSGHD.
. 
%%     ~~Describe \code{q} here~~
}
  \item{method}{ ( optional)  A string indicating the method used to pbtain output , if not specified MGHD is used. Alternative methods are:  MSGHD, cMSGHD, MCGHD. %%     ~~Describe \code{true.class} here~~
}
 \item{scale}{ ( optional) A  logical value indicating whether or not the data were scaled to obtain output, true by default.%%     ~~Describe \code{true.class} here~~
}}
\value{The contour plot
	}
\author{Cristina Tortora
Maintainer: Cristina Tortora <ctortora@mcmaster.ca>%%  ~~who you are~~
}
\examples{
##loading bankruptcy data
data(bankruptcy)

##model estimation
res=MCGHD(data=bankruptcy[,2:3],G=2,method="kmedoids",max.iter=30)
#result
contourpl(bankruptcy[,2:3],res,method="MCGHD")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

