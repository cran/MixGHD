\name{MSGHD}
\alias{MSGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture of multiple scaled generalized hyperbolic distributions (MSGHD).
%%  ~~function to do ... ~~
}
\description{Carries out model-based clustering using the mixture of multiple scaled  generalized hyperbolic distributions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{MSGHD(data=NULL,gpar0=NULL,G=2,max.iter=100,label=NULL,eps=1e-2,
	method="km",scale=TRUE,nr=10, modelSel="AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A n x p matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{x} here~~
}
  \item{gpar0}{ (optional) A list containing the initial parameters of the mixture model. See the 'Details' section. 
%%     ~~Describe \code{q} here~~
}
  \item{G}{ The range of values for the number of clusters.
%%     ~~Describe \code{G} here~~
}
  \item{max.iter}{(optional)  A numerical parameter giving  the maximum number of iterations each EM algorithm is allowed to use.
%%     ~~Describe \code{model} here~~
}
 
  \item{label}{ ( optional)  A n dimensional vector, if label[i]=k then observation i belongs to group k, If label[i]=0 then observation i has no known group, if NULL then the data has no known groups. %%     ~~Describe \code{true.class} here~~
}
 \item{eps}{ (optional) A number specifying the epsilon value for the convergence criteria used in the EM algorithms. For each algorithm, the criterion is based on the difference between the log-likelihood at an iteration and an asymptotic estimate of the log-likelihood at that iteration. This asymptotic estimate is based on the Aitken acceleration.
%%     ~~Describe \code{init} here~~
}
  \item{method}{ ( optional) A string indicating the initialization criterion, if not specified kmeans clustering is used. Alternative methods are: hierarchical "hierarchical", random "random", and model based "modelBased" clustering%%     ~~Describe \code{true.class} here~~
}
 \item{scale}{ ( optional) A  logical value indicating whether or not the data should be scaled, true by default.%%     ~~Describe \code{true.class} here~~
}
 \item{nr}{ ( optional) A  number indicating the number of starting value when random is used, 10 by default.%%     ~~Describe \code{true.class} here~~
}
\item{modelSel}{ ( optional) A string indicating the model selection criterion, if not specified AIC is used. Alternative methods are: BIC,ICL, and AIC3 %%     ~~Describe \code{true.class} here~~
}
}
\details{ The arguments gpar0, if specified, is a list structure containing at least one p dimensional vector  mu, alpha and phi, a pxp matrix gamma, and a px2 matrix cpl containing the vector omega and the vector lambda.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A S4 object of class \link{MixGHD} with slots:
	\item{index}{Value of the index used for model selection (AIC or ICL or BIC or AIC3) for each G,the index used is specified by the user, if not specified AIC is used.}
		\item{BIC}{Bayesian information criterion.}
	\item{ICL}{Integrated completed likelihood.}
	\item{AIC}{Akaike information criterion.}
	\item{AIC3}{Akaike information criterion 3.}
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
\references{C. Tortora, B.C. Franczak, R.P. Browne, and P.D. McNicholas (2019). A Mixture of Coalesced Generalized Hyperbolic Distributions. \emph{Journal of Classification} (to appear).
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora, Aisha ElSherbiny, Ryan P. Browne, Brian C. Franczak, and Paul D. McNicholas.
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{
 \code{\link{MGHD}}
}
\examples{
##loading banknote data
data(banknote)


##model estimation
model=MSGHD(banknote[,2:7],G=2,max.iter=30)

#result
table(banknote[,1],model@map)
summary(model)
plot(model)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Generalized hyperboilc distribution }

