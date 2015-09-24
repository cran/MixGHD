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
MGHD(data=NULL,gpar0=NULL,G=2,max.iter=100,label=NULL,eps=1e-2,
method="kmeans",scale=TRUE,nr=10, modelSel="AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A n x p matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{x} here~~
}
  \item{gpar0}{ (optional) A list containing the initial parameters of the mixture model. See the 'Details' section. 
%%     ~~Describe \code{q} here~~
}
  \item{G}{The range of values for the number of clusters.
%%     ~~Describe \code{G} here~~
}
  \item{max.iter}{(optional) A numerical parameter giving the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
   \item{label}{ ( optional)  A n dimensional vector, if label[i]=k then observation belongs to group k, If label[i]=0 then observation has no known group, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
 \item{eps}{ (optional) A number specifying the epsilon value for the convergence criteria used in the EM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration.
%%     ~~Describe \code{init} here~~
}
 \item{method}{ ( optional)  A string indicating the initialization criterion, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical", random "random", and model based "modelBased" clustering%%     ~~Describe \code{true.class} here~~
}
 \item{scale}{ ( optional) A  logical value indicating whether or not the data should be scaled, true by default.%%     ~~Describe \code{true.class} here~~
}
 \item{nr}{ ( optional) A  number indicating the number of starting value when random is used, 10 by default.%%     ~~Describe \code{true.class} here~~
}
\item{modelSel}{ ( optional) A string indicating the model selection criterion, if not specified AIC is used. Alternative methods are: BIC,ICL, and AIC3 %%     ~~Describe \code{true.class} here~~
}
}
\details{ The arguments gpar0, if specified, is a list structure containing at least one p dimensional vector  mu, and alpha, a pxp matrix sigma, and a 2 dimensional vector containing  omega and  lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with components
		\item{index}{Value of the index used for model selection (AIC or ICL or BIC or AIC3) for each G,the index used is specified by the user, if not specified AIC is used.}
	\item{model}{A list containing the following elements for the selected model}
	\item{BIC}{Bayesian information criterion value.}
	\item{ICL}{ICL index.}
		\item{AIC}{AIC index.}
	\item{AIC3}{AIC3 index.}
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
\references{ R.P. Browne, and  P.D. McNicholas (2015). A Mixture of  Generalized Hyperbolic Distributions. \emph{Canadian Journal of Statistics}
%% ~put references to the literature/web site here ~
}
\author{Ryan P. Browne and Cristina Tortora
Maintainer: Cristina Tortora <ctortora@mcmaster.ca>%%  ~~who you are~~
}
\examples{
##loading crabs data
data(crabs)

##model estimation
model=MGHD(data=crabs[,4:8],  G=2   )

#result
plot(model$model$loglik)
table(model$model$map, crabs[,2])

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
model=MGHD(data=bankruptcy[,2:3],  G=2, label=label )

#result
table(model$model$map,bankruptcy[,1])}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Classification }
\keyword{ Generalized hyperboilc distribution }% __ONLY ONE__ keyword per line

