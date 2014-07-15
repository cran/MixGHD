

MGHD <- function(data=NULL, gpar0=NULL, G=2, n=10, label =NULL   ) {
##Expexctation Maximization estimation of GHD
##data
## G n clusters
##n number of iterations
	data=scale(as.matrix(data))
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	if ( G < 1) stop('G is not a positive integer')
	if ( n < 1) stop('n is not a positive integer')
	
	if (is.null(gpar0)) gpar = igpar(data=data, g=G)
	else gpar  = gpar0
	
	loglik = numeric(n)
	for (i in 1:n) {
		gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label = label)	###parameter estimation	
		loglik[i] = llikGH(data, gpar) ##likelyhood
	}
	val = list(loglik= loglik, gpar=gpar, z=weightsGH(data=data, gpar= gpar), map=MAPGH(data=data, gpar= gpar, label=label) )
	return(val)
}