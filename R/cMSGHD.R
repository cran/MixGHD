
cMSGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL ,eps=1e-2, method="km" ,scale=TRUE) {
    data=as.matrix(data)
    if( scale==TRUE)
    {data=scale(data)}
        pcol=ncol(data)
        # if (nrow(data)<((G-1)+G*(4*pcol+pcol*(pcol-1)/2)))stop('G is too big, number of parameters > n')
    if (is.null(data)) stop('data is null')
    if (nrow(data) == 1) stop('nrow(data) is equal to 1')
    #  if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
    if (any(is.na(data))) stop('No NAs allowed.')
    if (is.null(G)) stop('G is NULL')
    if ( G < 1) stop('G is not a positive integer')
    if (  max.iter< 1) stop('max.iter is not a positive integer')
    n=max.iter
    if (is.null(gpar0)) gpar = rmgparMSr(g=G,data,method=method)
    else gpar  = gpar0
    
    loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label)
        loglik[i] = llikMS(data, gpar)
    }
   	while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
		i = i+1
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label)
        loglik[i] = llikMS(data, gpar)
    }
    if(i<n){loglik[i+1:max.iter]=loglik[i]}
    BIC=2*loglik[n]-log(nrow(data))*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    z=weightsMS(data=data, gpar= gpar)
    map=MAPMS(data=data, gpar= gpar, label=label)
    ICL=BIC+sum(log(apply(z,1,max)))
    val = list(loglik= loglik, gpar=gpar, z=z, map=map, BIC=BIC,ICL=ICL)
    return(val)
}