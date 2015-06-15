\name{ARI}
\alias{ARI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adjusted Rand Index.
%%  ~~function to do ... ~~
}
\description{Compares two classifications using the ddjusted Rand index (ARI).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ARI(x=NULL, y=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A n dimensional vector of class labels.
%%     ~~Describe \code{x} here~~
}
  \item{y}{  A n dimensional vector of class labels.
. 
%%     ~~Describe \code{q} here~~
}
}
\details{ The ARI has expected value 0 in case of random partition, it is equal to one in case of perfect agreement..
%%  ~~ If necessary, more details than the description above ~~
}
\value{The adjusted Rand index value
	}
\references{ L. Hubert and P. Arabie (1985) Comparing Partitions, \emph{Journal of the Classification} 2:193-218.
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora
Maintainer: Cristina Tortora <ctortora@mcmaster.ca>%%  ~~who you are~~
}
\examples{
##loading crabs data
data(banknote)

##model estimation
model=MGHD(data=banknote[,2:7],  G=2   )

#result
ARI(model$model$map, banknote[,1])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

