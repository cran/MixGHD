\name{MGHD}
\alias{MGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture of generalized hyperbolic distributions (MGHD).
%%  ~~function to do ... ~~
}
\description{Carries out model-based clustering and classification using the mixture of generalized hyperbolic distributions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MGHD(data=NULL, gpar0=NULL, G=2, n=10,  label=NULL )
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
  \item{n}{(optional) A numerical parameter giving the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
   \item{label}{ ( optional)  A n dimensional vector, if label[i]=k then observation belongs to group k, If label[i]=0 then observation has no known group, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
}
\details{ The arguments gpar0, if specified, is a list structure containing at least one p dimensional vector  mu, and alpha, a pxp matrix sigma, and a 1x2 matrix containing the vector omega and the vector lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with components
	\item{gpar}{A list of the model parameters. }
	\item{loglik}{The log-likelihood values.}
	\item{map}{A vector of integers indicating the maximum a posteriori classifications for the best model.}
	\item{z}{A matrix giving the raw values upon which map is based.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ R.P. Browne, and  P.D. McNicholas (2013). A Mixture of  Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1305.1036
%% ~put references to the literature/web site here ~
}
\author{Ryan P. Browne.
Maintainer: Cristina Tortora <ctortora@uoguelph.ca>%%  ~~who you are~~
}
\examples{
##loading crabs data
data(crabs)

##model estimation
model=MGHD(data=crabs[,4:8],  G=2, n=100   )

#result
plot(model$loglik)
table(model$map, crabs[,2])

## Classification
##loading bankruptcy data
data(bankruptcy)
#70% belong to the training set
 label=bankruptcy[,1]
#for a Classification porpuse the label cannot be 0
 label[1:33]=2
 a=round(runif(20)*65+1)
 label[a]=0
 
 
##model estimation
model=MGHD(data=bankruptcy[,2:3],  G=2, n=100  ,label=label )

#result
table(model$map,bankruptcy[,1])}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Classification }
\keyword{ Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

