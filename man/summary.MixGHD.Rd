\name{summary}
\alias{summary.MixGHD}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot objects of class MixGHD.
%%  ~~function to do ... ~~
}
\description{Methods for function \code{summary} aimed at summarizing the S4 classes included in the \code{\link{MixGHD}}-package
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A object of \code{\link{MixGHD-class}}.
%%     ~~Describe \code{x} here~~
}
  }\section{Methods}{
\describe{
\item{\code{signature(object = "MixGHD")}}{
S4 method for summaryzing objects of \code{\link{MixGHD-class}}.}
}
}
\author{Cristina Tortora.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{
 \code{\link{MixGHD}} \code{\link{MixGHD-class}},\code{\link{MGHD}},\code{\link{MCGHD}},\code{\link{MSGHD}},\code{\link{cMSGHD}},\code{\link{MGHFA}}
}
\examples{
##loading banknote data
data(bankruptcy)


##model estimation
model=MSGHD(bankruptcy[,2:3],G=2,max.iter=30)

#result

summary(model)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ methods}

