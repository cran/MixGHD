\name{MixGHD-class}
\alias{MixGHD-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Class MixGHD.
%%  ~~function to do ... ~~
}
\description{This class pertains to results of the application of function {\code{\link{MGHD}}},{\code{\link{MCGHD}}},{\code{\link{MSGHD}}},{\code{\link{cMSGHD}}}.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

%- maybe also 'usage' for other objects documented here.
\section{Slots}{
	\describe{
\item{Index}{Bayesian information criterion value for each combination of G and q.}
	\item{BIC}{Bayesian information criterion value.}
	\item{gpar}{A list of the model parameters. }
	\item{loglik}{The log-likelihood values.}
	\item{map}{A vector of integers indicating the maximum a posteriori classifications for the best model.}
	\item{z}{A matrix giving the raw values upon which map is based.}
	\item{method}{A string indicating the used method: MGHD, MGHFA, MSGHD, cMSGHD, MCGHD.}
	\item{data}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.}
	\item{par}{(only for MCGHD)A list of the model parameters in the rotated space.}
	}}
\details{ Plot the loglikhelyhood vale for each iteration of the EM algorithm. If p=2 it shows a contour plot. If 2<p<10 shows a splom of the data colored according to the cluster membership.
	}
\section{Methods}{
\describe{
\item{\code{signature(x = "MixGHD", y = "missing")}}{
S4 method for plotting objects of \code{\link{MixGHD-class}}.}
}
}
\author{Cristina Tortora, Aisha ElSherbiny, Ryan P. Browne, Brian C. Franczak, and Paul D. McNicholas.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{
 \code{\link{MixGHD-class}},\code{\link{MGHD}},\code{\link{MCGHD}},\code{\link{MSGHD}},\code{\link{cMSGHD}},\code{\link{MGHFA}}
}
\examples{
##loading banknote data
data(bankruptcy)


##model estimation
model=MSGHD(bankruptcy[,2:3],G=2,max.iter=30)

#result
summary(model)
plot(model)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Plot}

