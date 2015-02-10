

MGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL , eps=1e-2, method="kmeans" ,scale=TRUE ) {
##Expexctation Maximization estimation of GHD
##data
## G n clusters
##n number of iterations
  data=as.matrix(data)
if( scale==TRUE){
	data=scale(as.matrix(data))}
    pcol=ncol(data)
    #if (nrow(data)<((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2)))stop('G is too big, number of parameters > n')
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	if (any(is.na(data))) stop('No NAs allowed.')
	if (is.null(G)) stop('G is NULL')
	if ( G < 1) stop('G is not a positive integer')
	if ( max.iter < 1) stop('max.iter is not a positive integer')
	
	if (is.null(gpar0)) gpar = igpar(data=data, g=G, method=method)
	else gpar  = gpar0
    
    n=max.iter
		loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label = label)	###parameter estimation
		loglik[i] = llikGH(data, gpar)}
        
        while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
            i = i+1
		gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label = label)	###parameter estimation

		loglik[i] = llikGH(data, gpar) ##likelyhood
	}

    BIC=2*loglik[n]-log(nrow(data))*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
	val = list(loglik= loglik, gpar=gpar, z=weightsGH(data=data, gpar= gpar), map=MAPGH(data=data, gpar= gpar, label=label),BIC=BIC )
	return(val)
}
