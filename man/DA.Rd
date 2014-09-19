\name{DA}
\alias{DA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Discriminant analysis using the mixture of  generalized hyperbolic distributions.
%%  ~~function to do ... ~~
}
\description{Carries out model-based discriminant analysis using 4 different models: the mixture of multiple scaled  generalized hyperbolic distributions (MGHD), the mixture of generalized hyperbolic factor analyzers (MGHFA), the mixture of multiple scaled  generalized hyperbolic distributions (MSGHD) and the mixture of coaelesed generalized hyperbolic distributions (MCGHD).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DA(training,traininglabel,test,testlabel, method="GHD",starting="kmeans",max.iter=10,q=2)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training}{ A n1 x p matrix or data frame such that rows correspond to observations and columns correspond to variables of the training data set.
%%     ~~Describe \code{x} here~~
}
\item{traininglabel}{ A n1 dimensional vector of membership for the units of the training set. If label[i]=k then observation belongs to group k. %%     ~~Describe \code{true.class} here~~
}
  \item{test}{ A n2 x p matrix or data frame such that rows correspond to observations and columns correspond to variables of the test data set.
%%     ~~Describe \code{x} here~~
}
\item{testlabel}{ A n2 dimensional vector of membership for the units of the test set. If label[i]=k then observation belongs to group k. %%     ~~Describe \code{true.class} here~~
}
  \item{method}{ ( optional)  A string indicatid the method to be used form discriminant analysis , if not specified GHD is used. Alternative methods are: MGHFA, MSGHD, MCGHD. %%     ~~Describe \code{true.class} here~~
}
  
  \item{starting}{ ( optional)  A string indicatid the initializtion criteris, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical" and model based "modelBased" %%     ~~Describe \code{true.class} here~~
}
\item{max.iter}{(optional)  A numerical parameter giving  the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
  \item{q}{(optional) used only if MGHFA method is selected. A numerical parameter giving the number of factors.
%%     ~~Describe \code{G} here~~
}
}
\value{A list with components
	\item{model}{A list with the model parameters.}
	\item{testMembership}{A vector of integers indicating the membership of the units in the test set}
	\item{ARI }{A value indicating the adjusted rand index for the test set.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{R.P. Browne, and  P.D. McNicholas (2013). A Mixture of  Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1305.1036\cr
C. Tortora, B.C. Franczak, R.P. Browne, and  P.D. McNicholas (2014). A Mixture of Coalesced Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1403.2332 \cr
C. Tortora,  P.D. McNicholas, and R.P. Browne (2014). A Mixture of Generalized Hyperbolic Factor Analyzers. Arxiv preprint arXiv:1311.6530%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora
	Maintainer: Cristina Tortora <ctortora@mcmaster.ca>%%  ~~who you are~~
}
\seealso{
 \code{\link{MGHD}} \code{\link{MGHFA}} \code{\link{MSGHD}} \code{\link{MCGHD}}
}
\examples{
##loading crabs data
data(crabs)


##divide the data in training set and test set
lab=rbind(matrix(1,100,1),matrix(2,100,1))
train=crabs[c(1:90,110:200),4:8]
trainl=lab[c(1:90,110:200),1]
test=crabs[91:109,4:8]
testl=lab[91:109,1]


##model estimation
 model=DA(train,trainl,test,testl,method="MSGHD",max.iter=40,starting="hierarchical")

#result
model$ARItest}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Generalized hyperboilc distribution }

