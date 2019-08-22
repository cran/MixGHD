\name{DA}
\alias{DA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Discriminant analysis using the mixture of  generalized hyperbolic distributions.}
%%  ~~function to do ... ~~

\description{Carries out model-based discriminant analysis using 5 different models: the mixture of multiple scaled  generalized hyperbolic distributions (MGHD), the mixture of generalized hyperbolic factor analyzers (MGHFA), the mixture of multiple scaled  generalized hyperbolic distributions (MSGHD),the mixture of convex multiple scaled  generalized hyperbolic distributions (cMSGHD) and the mixture of coaelesed generalized hyperbolic distributions (MCGHD).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{DA(train,trainL,test,testL,method="MGHD",starting="km",max.iter=100,
	eps=1e-2,q=2,scale=TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train}{ A n1 x p matrix or data frame such that rows correspond to observations and columns correspond to variables of the training data set.
%%     ~~Describe \code{x} here~~
}
\item{trainL}{ A n1 dimensional vector of membership for the units of the training set. If trainL[i]=k then observation belongs to group k. %%     ~~Describe \code{true.class} here~~
}
  \item{test}{ A n2 x p matrix or data frame such that rows correspond to observations and columns correspond to variables of the test data set.
%%     ~~Describe \code{x} here~~
}
\item{testL}{ A n2 dimensional vector of membership for the units of the test set. If testL[i]=k then observation belongs to group k. %%     ~~Describe \code{true.class} here~~
}
  \item{method}{ ( optional)  A string indicating the method to be used form discriminant analysis , if not specified MGHD is used. Alternative methods are: MGHFA, MSGHD, cMSGHD, MCGHD. %%     ~~Describe \code{true.class} here~~
}
   \item{starting}{ ( optional) A string indicating the initialization criterion, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical",random "random", kmedoids "kmedoids", and model based "modelBased" %%     ~~Describe \code{true.class} here~~
}
\item{max.iter}{(optional)  A numerical parameter giving  the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
 \item{eps}{ (optional) A number specifying the epsilon value for the convergence criteria used in the EM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration.
%%     ~~Describe \code{init} here~~
}
  \item{q}{(optional) used only if MGHFA method is selected. A numerical parameter giving the number of factors.
%%     ~~Describe \code{G} here~~
}
  \item{scale}{ ( optional) A  logical value indicating whether or not the data should be scaled, true by default.%%     ~~Describe \code{true.class} here~~
}
}
\value{A list with components
	\item{model}{An S4 object of class \code{\link{MixGHD}} with the model parameters.}
	\item{testMembership}{A vector of integers indicating the membership of the units in the test set}
	\item{ARItest }{A value indicating the adjusted rand index for the test set.}
	\item{ARItrain }{A value indicating the adjusted rand index for the train set.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{R.P. Browne, and  P.D. McNicholas (2015). A Mixture of  Generalized Hyperbolic Distributions. \emph{Canadian Journal of Statistics}, 43.2 176-198.\cr
C. Tortora, B.C. Franczak, R.P. Browne, and P.D. McNicholas (2019). A Mixture of Coalesced Generalized Hyperbolic Distributions. \emph{Journal of Classification} (to appear). \cr
C.Tortora, P.D. McNicholas, and R.P. Browne (2016). Mixtures of Generalized Hyperbolic Factor Analyzers. 
\emph{Advanced in data analysis and classification} 10(4) p.423-440.%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora, Aisha ElSherbiny, Ryan P. Browne, Brian C. Franczak, and Paul D. McNicholas.
	Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{ \code{"\link[=MixGHD-class]{MixGHD}"} 
 \code{\link{MGHD}} \code{\link{MGHFA}} \code{\link{MSGHD}} \code{\link{cMSGHD}} \code{\link{MCGHD}} \code{\link{ARI}}  \code{\link{MixGHD-class}} \code{\link{MixGHD}}
}
\examples{
##loading banknote data
data(banknote)
banknote[,1]=as.numeric(factor(banknote[,1]))


##divide the data in training set and test set
train=banknote[c(1:74,126:200),]
test=banknote[75:125,]

##model estimation
 model=DA(train[,2:7],train[,1],test[,2:7],test[,1],method="MGHD",max.iter=20)

#result
model$ARItest}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Generalized hyperboilc distribution }

