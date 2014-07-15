\name{sonar}
\alias{sonar}
\docType{data}

\title{Sonar data}

\description{ The data report the patterns obtained by bouncing sonar signals at various angles and under various conditions. There are 208 patterns in all, 111 obtained by bouncing sonar signals off a metal cylinder and 97 obtained by bouncing signals off rocks. Each pattern is a set of 60 numbers (variables) taking values between 0 and 1.}


\usage{data(sonar)}

\format{A data frame with 208 observations and 61 columns. The first 60 columns contain the variables. The 61st column gives the material: \code{1} rock), \code{2} metal.}
\source{UCI machine learning repository}
\references{R.P. Gorman and T. J. Sejnowski (1988) Analysis of hidden units in a layered network trained to classify sonar targets.
 \emph{Neural Networks} \bold{1}: 75-89
}
% \examples{}
\keyword{datasets}
