mainMGHFA<-function(data=NULL, gpar0, G, n, label  , eps, method ,q,nr=nr ) {
    pcol=ncol(data)
    if(!is.null(label)){
      lc=apply(data[label==1,],2,mean)
    #  if(min(label)==0&max(label)==G){
        for(i in 2:G){
          lc=rbind(lc,apply(data[label==i,],2,mean))
        }#}
   
      
      z = combinewk(weights=matrix(1/G,nrow=nrow(data),ncol=G), label=label)
      if (is.null(gpar0)) gpar  = rgpar(data=data, g=G, w=z,l=lc)
      else gpar  = gpar0
          }
    else{
        if (is.null(gpar0)) gpar = igparM(data=data, g=G,q=q,method=method,nr=nr)
        else gpar  = gpar0}

        loglik = numeric(n)
        for (i in 1:3) {
                gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)     ###parameter estimation
                loglik[i] = llikFA(data, gpar) ##likelyhood
        }
        while ( ( getall(loglik[1:i]) > eps) & (i < (n) ) )  {
          i = i+1
          gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)	###parameter estimation
          loglik[i] = llikFA(data, gpar) ##likelyhood
        }


#    if(i<n){loglik=loglik[-(i+1:n)]}
     if(i<n){loglik=loglik[-(i+1:n)]}
     BIC=2*loglik[i]-log(nrow(data))*((G-1)+G*(3*pcol+2+pcol*q-q*(q-1)/2))
     val = list(loglik= loglik, gpar=gpar, z=weightsFA(data=data, gpar= gpar), map=MAPFA(data=data, gpar= gpar, label=label) , BIC=BIC)
     return(val)
}




MGHFA<- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL  ,q=2,eps=1e-2, method="kmeans",scale=TRUE ,nr=10) {
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
    #	if ( G < 1) stop('G is not a positive integer')
	if ( max.iter < 1) stop('max.iter is not a positive integer')
	#if ( q < 1) stop('n is not a positive integer')
    bico=-Inf
    t=length(G)
    tq=length(q)
    BIC=matrix(NA,t,tq)
    cont=0
    if(length(G)==1&length(q)==1){
        mo=(mainMGHFA(data=data, gpar0=gpar0, G=G,q=q, n=max.iter, eps=eps,  label=label,method= method,nr=nr))#,silent = TRUE
        model=mo
        val=MixGHD(BIC=model$BIC, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHFA",data=as.data.frame(data))
        
    }
    else{
    for(b2 in 1:tq){
	for(b in 1:t){
        ct=1
        while(ct<4 || is.list(mo)==FALSE){
            mo=try(mainMGHFA(data=data, gpar0=gpar0, G=G[b],q=q[b2], n=max.iter, eps=eps,  label=label,method= method),silent = TRUE)
        ct=ct+1}
        cont=cont+1
        if(is.list(mo)){
            bicn=mo$BIC
            BIC[b,b2]=bicn}
        else{bicn=-Inf
            BIC[b,b2]=NA}
        if(bicn>bico){
            bico=bicn
            sg=G[b]
            sq=q[b2]
            model=mo
        }
    }}
      val=MixGHD(Index=BIC,BIC=model$BIC, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHFA",data=as.data.frame(data),scale=scale)
    

        cat("The best model (BIC) for the range of factors and components used is  G = ", sg,", and q=", sq ,".\nThe BIC for this model is ", bico,".",sep="")}

    return(val)


}



