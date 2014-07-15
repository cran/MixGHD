\name{MixGHD-package}
\alias{MixGHD-package}
\alias{MixGHD}
\docType{package}
\title{Model based clustering and classification using the mixture of generalized hyperbolic distributions

}
\description{Carries out model-based clustering using three different models. The models are all based on the  generalized hyperbolic distribution.The first model MGHD is the classical mixture of generalized hyperbolic distributions. The MGHFA is the mixture of generalized hyperbolic factor analyzers for high dimensional data sets. The MCGHD, mixture of  coalesced generalized hyperbolic distributions is a new more flexible model.
}
\details{
\tabular{ll}{
Package: \tab MixGHD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-19\cr
License: \tab GPL (>=2) \cr
}
This package contains the function MGHD for model based clustering and classification using the mixture of generalized hyperbolic distributions;
 the function MGHFA for model based clustering and classification using the mixture of generalized hyperbolic factor analyzers;
the function MCGHD for model based clustering  using the mixture of coalesced generalized hyperbolic distributions, and some real data sets.
}
\author{
Cristina Tortora, Ryan P. Browne, Brian C. Franczak and Paul D. McNicholas.
Maintainer: Cristina Tortora <ctortora@uoguelph.ca>%%  ~~who you are~~
}
\references{
R.P. Browne, and  P.D. McNicholas (2013). A Mixture of  Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1305.1036\cr
C. Tortora, B.C. Franczak, R.P. Browne, and  P.D. McNicholas (2014). A Mixture of Coalesced Generalized Hyperbolic Distributions. Arxiv preprint arXiv:1403.2332 \cr
C. Tortora,  P.D. McNicholas, and R.P. Browne (2014). A Mixture of Generalized Hyperbolic Factor Analyzers. Arxiv preprint arXiv:1311.6530
}
\keyword{ Clustering }
\keyword{ Classification }
\keyword{ Expectation Maximization algorithm }
\keyword{ Generalized hyperboilc distribution }

\examples{
###model based clustering  using the mixture of generalized hyperbolic distributions
data(crabs)
modelGH=MGHD(data=crabs[,4:8],  G=2, n=50)
table(modelGH$map,crabs[,2])

##model based classification using the mixture of generalized hyperbolic factor analyzers
## Classification
#70% belong to the training set
data(sonar)
 label=sonar[,61]
 a=round(runif(62)*207+1)
 label[a]=0
 
model=MGHFA(data=sonar[,1:60],  G=2, n=30  ,q=2,label=label )
table(model$map,sonar[,61])


##model based clustering  using the mixture of generalized hyperbolic factor analyzers
data(bankruptcy)
 modelCGH=MCGHD(bankruptcy[,2:3],G=2)
 table(modelCGH$map,bankruptcy[,1])

}
