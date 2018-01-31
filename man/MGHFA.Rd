\name{MGHFA}
\alias{MGHFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture of generalized hyperbolic factor analyzers (MGHFA).
%%  ~~function to do ... ~~
}
\description{Carries out model-based clustering and classification using the mixture of generalized hyperbolic factor analyzers.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MGHFA(data=NULL, gpar0=NULL, G=2, max.iter=100, 
label =NULL  ,q=2,eps=1e-2 , method="kmeans", scale=TRUE ,nr=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{x} here~~
}
  \item{gpar0}{ (optional) A list containing the initial parameters of the mixture model. See the 'Details' section. 
%%     ~~Describe \code{q} here~~
}
  \item{G}{ The range of values for the number of clusters.
%%     ~~Describe \code{G} here~~
}
  \item{max.iter}{(optional) A numerical parameter giving the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
  \item{label}{ ( optional)  A n dimensional vector, if label[i]=k then observation i belongs to group k, If label[i]=0 then observation i has no known group, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
 \item{q}{ The range of values for the number of  factors.
%%     ~~Describe \code{G} here~~
}
 \item{eps}{ (optional) A number specifying the epsilon value for the convergence criteria used in the EM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration.
%%     ~~Describe \code{init} here~~
}
 \item{method}{ ( optional)  A string indicating the initialization criterion, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical" and model based "modelBased" clustering %%     ~~Describe \code{true.class} here~~
}
 \item{scale}{ ( optional) A  logical value indicating whether or not the data should be scaled, true by default.%%     ~~Describe \code{true.class} here~~
}
 \item{nr}{ ( optional) A  number indicating the number of starting value when random is used, 10 by default.%%     ~~Describe \code{true.class} here~~
}
}
\details{ The arguments gpar0, if specified, is a list structure containing at least one p dimensional vector  mu, alpha and phi, a pxp matrix gamma, a 2 dimensional vector cpl containing omega and lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A S4 object of class \link{MixGHD}  with slots:
	\item{Index}{Bayesian information criterion value for each combination of G and q.}
	\item{BIC}{Bayesian information criterion value.}
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
\references{C.Tortora, P.D. McNicholas, and R.P. Browne. Mixtures of Generalized Hyperbolic Factor Analyzers. 
\emph{Advanced in data analysis and classification} 10(4) p.423-440, 2016.
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora, Ryan P. Browne, and Paul D. McNicholas.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\examples{
## Classification
#70% belong to the training set
data(sonar)
 label=sonar[,61]
 set.seed(4)
 a=round(runif(62)*207+1)
 label[a]=0
 
 
##model estimation
model=MGHFA(data=sonar[,1:60],  G=2, max.iter=25  ,q=2,label=label)

#result
table(model@map,sonar[,61])
summary(model)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Classification }
\keyword{ Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

