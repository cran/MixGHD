\name{Selection}
\alias{Selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of parameter using the mixture of  generalized hyperbolic distributions.
%%  ~~function to do ... ~~
}
\description{Replicate a routine n times and gives as result the best solution, the available routine are:
the mixture of multiple scaled  generalized hyperbolic distributions (MGHD), the mixture of generalized hyperbolic factor analyzers (MGHFA), the mixture of multiple scaled  generalized hyperbolic distributions (MSGHD), the mixture of convex multiple scaled  generalized hyperbolic distributions (cMSGHD) and the mixture of coaelesed generalized hyperbolic distributions (MCGHD).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{Selection(data,label, G=2, niter=50,method="GHD",starting="kmeans",
	max.iter=10,eps=1e-2,labelcl=NULL,q=2,scale=TRUE,criterion="ARI")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A n x p matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{x} here~~
}
\item{label}{  A n dimensional vector of known membership, if label[i]=k unit i belongs to cluster k. It is used to select the solution.   %%     ~~Describe \code{true.class} here~~
}
   \item{G}{ A numerical parameter giving the number of clusters.
%%     ~~Describe \code{G} here~~
}
   \item{niter}{ A numerical parameter giving the number of replications.
%%     ~~Describe \code{G} here~~
}
  \item{method}{ ( optional)  A string indicating the method to be used, if not specified GHD is used. Alternative methods are: MGHFA, MSGHD, cMSGHD, MCGHD. %%     ~~Describe \code{true.class} here~~
}
  
  \item{starting}{ ( optional)  A string indicating the initialization criterion, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical",random "random", and model based "modelBased" %%     ~~Describe \code{true.class} here~~
}
\item{max.iter}{(optional)  A numerical parameter giving  the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
 \item{eps}{ (optional) A number specifying the epsilon value for the convergence criteria used in the EM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration.
%%     ~~Describe \code{init} here~~
}
\item{labelcl}{ (optional) A n dimensional vector of membership, if labelcl[i]=k then observation belongs to group k, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
  \item{q}{(optional) used only if MGHFA method is selected. A numerical parameter giving the number of factors.
%%     ~~Describe \code{G} here~~
}
 \item{scale}{ ( optional) A  logical value indicating whether or not the data should be scaled, true by default.%%     ~~Describe \code{true.class} here~~
}
 \item{criterion}{ ( optional) A string indicating the selection criterion, if not specified the adjusted rand index (ARI)  is used. The alternative criteria are the BIC and the ICL.%%     ~~Describe \code{true.class} here~~
}
}
\details{ Select the best model according to the adjusted rand index or the BIC
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with components
	\item{model}{A list with the model parameters.}
	\item{ARI }{A value indicating the adjusted rand index.}
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
 \code{\link{MGHD}} \code{\link{MGHFA}} \code{\link{MSGHD}} \code{\link{cMSGHD}} \code{\link{MCGHD}}
}
\examples{
##loading bankruptcy data
#data(bankruptcy)

##model estimation
#mod=Selection(bankruptcy[,2:3],bankruptcy[,1],method="MCGHD")

##result
#mod$ARI}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Generalized hyperboilc distribution }

