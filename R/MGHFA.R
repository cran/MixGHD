
MGHFA<- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL  ,q=2,eps=1e-2, method="kmeans",scale=TRUE ) {
##Expexctation Maximization estimation of GHD
##data
## G n clusters
##n number of iterations
  data=as.matrix(data)
if( scale==TRUE)
{data=scale(data)}
	pcol=ncol(data)
    #  if (nrow(data)<((G-1)+G*(3*pcol+2+pcol*q-q*(q-1)/2)))stop('G is too big, number of parameters > n')
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function  only works with multivariate data p > 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	if ( G < 1) stop('G is not a positive integer')
	if ( max.iter < 1) stop('max.iter is not a positive integer')
	if ( q < 1) stop('n is not a positive integer')	
	n=max.iter
	if (is.null(gpar0)) gpar = igparM(data=data, g=G,q=q,method=method)
	else gpar  = gpar0
	loglik = numeric(n)
	for (i in 1:3) {
		gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)	###parameter estimation	
		loglik[i] = llikFA(data, gpar) ##likelyhood
	}
	while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
		i = i+1
		gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)	###parameter estimation	
		loglik[i] = llikFA(data, gpar) ##likelyhood
	}
    if(i<n){loglik[i+1:max.iter]=loglik[i]}
	BIC=2*loglik[n]-log(nrow(data))*((G-1)+G*(3*pcol+2+pcol*q-q*(q-1)/2))
	val = list(loglik= loglik, gpar=gpar, z=weightsFA(data=data, gpar= gpar), map=MAPFA(data=data, gpar= gpar, label=label) , BIC=BIC)
	return(val)
}



