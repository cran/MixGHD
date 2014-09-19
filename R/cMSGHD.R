
cMSGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL , method="kmeans" ) {
    data=as.matrix(data)
    data=scale(data)
    if (is.null(data)) stop('data is null')
    if (nrow(data) == 1) stop('nrow(data) is equal to 1')
    if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
    if (any(is.na(data))) stop('No NAs allowed.')
    if (is.null(G)) stop('G is NULL')
    if ( G < 1) stop('G is not a positive integer')
    if (  max.iter< 1) stop('max.iter is not a positive integer')
    n=max.iter
    if (is.null(gpar0)) gpar = rgparMSr(g=G,p=ncol(data),data,method=method)
    else gpar  = gpar0
    
    loglik = numeric(n)
    for (i in 1:n) {
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label)
        loglik[i] = llikMS(data, gpar)
    }
    pcol=ncol(data)
    BIC=2*loglik[n]-log(nrow(data))*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    val = list(loglik= loglik, gpar=gpar, z=weightsMS(data=data, gpar= gpar), map=MAPMS(data=data, gpar= gpar, label=label), BIC=BIC)
    return(val)
}