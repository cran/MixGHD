
MCGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, epsilon=1e-2,  label=NULL ) {
	data=as.matrix(data)
	data=scale(data)
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	if ( G < 1) stop('G is not a positive integer')
	if (  max.iter< 1) stop('n is not a positive integer')
	
	
	if (is.null(gpar0)) gpar  = rmgpar(g=G,p=ncol(data),data=data)
	else gpar = gpar0
	
	loglik = numeric(max.iter)
	for (i in 1:3) {
		gpar = EMgrstep(data=data, gpar=gpar, v=1, label = label)	
		loglik[i] = llik(data, gpar)
	}
	
	while ( ( getall(loglik[1:i]) > epsilon) & (i < (max.iter) ) )  {
		i = i+1
		gpar = EMgrstep(data=data, gpar=gpar, v=1, label = label)		
		loglik[i] = llik(data, gpar)
		
	}
	pcol=ncol(data)
	BIC=2*loglik[max.iter]-log(nrow(data))*(2*(G-1)+G*(2*pcol+0.5*pcol*(pcol-1))+G*2*pcol+G*2)
#BIC=2*loglik[max.iter]-log(nrow(data))*((G-1)+G*(3*pcol+0.5*pcol*(pcol-1)))
	par=partrue(gpar,G)
	val = list(loglik= loglik[1:i], gpar=gpar,par=par, z=weights(data=data, gpar= gpar), map=MAP(data=data, gpar= gpar, label=label),BIC=BIC )
	return(val)
}

