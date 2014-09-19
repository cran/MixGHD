\name{rMSGHD}
\alias{cMSGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convex mixture of multiple scaled generalized hyperbolic distributions (MSGHD).
%%  ~~function to do ... ~~
}
\description{Carries out model-based clustering using the convex mixture of multiple scaled  generalized hyperbolic distributions. The cMSGHD only allows conves level sets.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cMSGHD(data=NULL, gpar0=NULL, G=2, max.iter=100, label=NULL, method="kmeans" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A n x p matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{x} here~~
}
  \item{gpar0}{ (optional) A list containing the initial parameters of the mixture model. See the 'Details' section. 
%%     ~~Describe \code{q} here~~
}
  \item{G}{ A numerical parameter giving the number of clusters.
%%     ~~Describe \code{G} here~~
}
  \item{max.iter}{(optional)  A numerical parameter giving  the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
 
  \item{label}{ ( optional)  A n dimensional vector, if label[i]=k then observation belongs to group k, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
  \item{method}{ ( optional) AA string indicating the initialization criteria, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical" and model based "modelBased" %%     ~~Describe \code{true.class} here~~
}
}
\details{ The arguments gpar0, if specified, is a list structure containing at least one p dimensional vector  mu, alpha and phi, a pxp matrix gamma, and a px2 matrix containing the vector omega0 and the vector lambda0.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with components
	\item{BIC}{Bayesian information criterion value.}
	\item{gpar}{A list of the model parameters }
	\item{loglik}{The log-likelihood values.}
	\item{map}{A vector of integers indicating the maximum a posteriori classifications for the best model.}
	\item{z}{A matrix giving the raw values upon which map is based.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{C. Tortora, B.C. Franczak, R.P. Browne, and  P.D. McNicholas (2014). A Mixture of Coalesced Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1403.2332
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora, Ryan P. Browne, and Paul D. McNicholas.
Maintainer: Cristina Tortora <ctortora@mcmaster.ca>%%  ~~who you are~~
}
\seealso{
 \code{\link{MGHD}} \code{\link{MSGHD}}
}
\examples{
##loading banknote data
data(banknote)

##model estimation
model=cMSGHD(banknote[,2:7],G=2,max.iter=30)

#result
table(banknote[,1],model$map)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Generalized hyperboilc distribution }

