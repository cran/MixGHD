\name{MixGHD class}
\alias{MixGHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Class "MixGHD"
%%  ~~function to do ... ~~
}
\description{This class pertains to results of the application of function \link{MGHD}, \link{MSGHD}, \link{cMSGHD}, \link{MCGHD}, and \link{MGHFA}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\section{Objects from the Class}{
Objects can be created  as a result to a call to \link{MGHD}, \link{MSGHD}, \link{cMSGHD}, \link{MCGHD}, and \link{MGHFA}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
  	\item{\code{index}}{Value of the index used for model selection (AIC or ICL or BIC or AIC3) for each G,the index used is specified by the user, if not specified AIC is used.}
		\item{BIC}{Bayesian information criterion value.}
	\item{\code{ICL}}{ICL index.}
	\item{\code{AIC}}{AIC index.}
	\item{\code{AIC3}}{AIC3 index.}
	\item{\code{gpar}}{A list of the model parameters (in the rotated space for MCGHD). }
	\item{\code{loglik}}{The log-likelihood values.}
	\item{\code{map}}{A vector of integers indicating the maximum a posteriori classifications for the best model.}
	\item{\code{par}}{Only for MCGHD. A list of the model parameters.}
	\item{\code{z}}{A matrix giving the raw values upon which map is based.}

   }}
\section{Methods}{
  \describe{
        	\item{plot}{\code{signature(x = "MixGHD")} 
	Provides plots of \code{\link{MixGHD-class}} by plotting
    the following elements: 	

	\itemize{
		\item {the value of the log likelihood for each iteration. }
		\item {Scatterplot of the data of all the possible couples of coordinates coloured according to the cluster. Only for less than 10 variables.}
		\item {If the number of variables is two: scatterplot and contour plot of the data coloured according to the cluster}}}
		

    \item{summary}{\code{summary(x = "MixGHD")}.

	Provides a summary of \code{\link{MixGHD-class}} objects by printing
    the following elements: 	
		\itemize{
		\item {The number   components used for the model}
		\item {BIC;}
		\item {AIC;}
		\item {AIC3;}
		\item {ICL;}
		\item {A table with the number of element in each cluster.}
		} }
	 
}}

\author{Cristina Tortora
Maintainer: Cristina Tortora <cristina.tortora@sjsu.edu>%%  ~~who you are~~
}
\seealso{
 \code{\link{MixGHD-class}}}
\examples{
##loading bankruptcy data
data(bankruptcy)

##model estimation
#res=MCGHD(data=bankruptcy[,2:3],G=2,method="kmedoids",max.iter=30)
#result
#plot(res)
#summary(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

