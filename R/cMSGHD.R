     run_EMstepMSr <- function(data = NULL, gpar = NULL, loglik=NULL, maxit=NULL, N= NULL, p=NULL, G=NULL, eps=NULL, label=NULL){
     counter = 0
     mu = matrix(0,nrow=G,ncol=p)
     phi= matrix(0,nrow=G,ncol=p)
     alpha= matrix(0,nrow=G,ncol=p)
     cplr = matrix(0,nrow=G,ncol=p*2)
     gamr = matrix(0,nrow=G,ncol=p*p)
     for(k in 1:G){
     par = gpar[[k]]
     mu[k,] = par$mu
     phi[k,] = par$phi
     alpha[k,]= par$alpha
     cplr[k,] = c(t(par$cpl))
     gamr[k,] = c(t(par$gam))
     }
     pi = gpar$pi
     v = 1
     mu1=t(mu)
     phi1=t(phi)
     alpha1=t(alpha)
     gam1=t(gamr)
     cpl1=t(cplr)
     if(is.null(label) ==T) label=rep(0,N)
     EMSr <- .C("EMstepSr", as.double(mu1), as.double(phi1), as.double(alpha1),
             as.double(cpl1), as.double(gam1), as.double(data), as.double(loglik), as.integer(maxit), as.integer(N), as.integer(p), as.integer(G),
as.double(eps), as.integer(label), as.double(pi), as.integer(v),
         as.integer(counter), PACKAGE="MixGHD")
     loglik = EMSr[[7]]
     counter = EMSr[[16]]
     mu = matrix(EMSr[[1]],nrow =G, ncol=p, byrow=TRUE)
     phi= matrix(EMSr[[2]], nrow=G, ncol=p, byrow=TRUE)
     alpha= matrix(EMSr[[3]], nrow=G, ncol=p, byrow=TRUE)
     cpl = array(EMSr[[4]],dim=c(p,2,G))
     gam = array((EMSr[[5]]),dim=c(p,p,G))
     gpar = list()
     for(k in 1:G){
        gpar[[k]] = list()
        gpar[[k]]$mu = mu[k,]
        gpar[[k]]$phi= phi[k,]
        gpar[[k]]$alpha = alpha[k,]
        gpar[[k]]$cpl = matrix(cpl[,,k], nrow=p, ncol=2, byrow=TRUE)
        gpar[[k]]$gam = matrix(gam[,,k], nrow=p, ncol=p, byrow=TRUE)
     }
     gpar$pi = EMSr[[14]]
     val = list(loglik,gpar,counter)
     return(val)
}

maincMSGHD<-function(data=NULL, gpar0=NULL, G, n, label  ,eps, method,nr=NULL ) {
    pcol=ncol(data)
    if(!is.null(label)){
      lc=apply(data[label==1,],2,mean)
   #   if(min(label)==0&max(label)==G){
        for(i in 2:G){
          lc=rbind(lc,apply(data[label==i,],2,mean))
        }#}
  
      
      z = combinewk(weights=matrix(1/G,nrow=nrow(data),ncol=G), label=label)
      if (is.null(gpar0)) gpar  = rgparMSr(data=data, g=G, w=z,l=lc)
      else gpar  = gpar0
      
    }
  
    else{
        if (is.null(gpar0)) gpar = rmgparMSr(g=G,data,method=method,nr=nr)
        else gpar  = gpar0}


    loglik = numeric(n)
    for (i in 1:3) {
        gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llikMS(data, gpar)
    }
       N = nrow(data)
       p = ncol(data)
       maxit = n
       temp <- run_EMstepMSr(data, gpar, loglik, maxit, N, p, G, eps, label)
       i = temp[[3]]


    if(i<n){temp[[1]]=temp[[1]][-(i+1:n)]}
    BIC=2*temp[[1]][i]-log(nrow(data))*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    z=weightsMS(data=data, gpar= temp[[2]])
    map=MAPMS(data=data, gpar= temp[[2]], label=label)
    ICL=BIC+2*sum(log(apply(z,1,max)))
    AIC=2*temp[[1]][i]-2*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    AIC3=2*temp[[1]][i]-3*((G-1)+G*(4*pcol+pcol*(pcol-1)/2))
    val = list(loglik= temp[[1]], gpar=temp[[2]], z=z, map=map, BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3)
    return(val)

}


cMSGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL ,eps=1e-2, method="km" ,scale=TRUE,nr=10, modelSel="AIC") {
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
    #if ( G < 1) stop('G is not a positive integer')
    if (  max.iter< 1) stop('max.iter is not a positive integer')

   

if(modelSel=="BIC"){


bico=-Inf
t=length(G)
BIC=matrix(NA,t,1)
cont=0
for(b in 1:t){
    mo=try(maincMSGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label, method=method,nr=nr),silent = TRUE)
    cont=cont+1
    if(is.list(mo)){
        bicn=mo$BIC
        BIC[cont]=bicn}
    else{bicn=-Inf
        BIC[cont]=NA}
    if(bicn>bico){
        bico=bicn
        sg=G[b]
        model=mo
    }
}
# val=list(BIC=BIC,model=model)
val=MixGHD(Index=BIC,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z, method="cMSGHD",data=as.data.frame(data),scale=scale)
cat("The best model (BIC) for the range of  components used is  G = ", sg,".\nThe BIC for this model is ", bico,".",sep="")
    return(val)}



else if(modelSel=="ICL"){
    bico=-Inf
    t=length(G)
    ICL=matrix(NA,t,1)
    cont=0
    for(b in 1:t){
        mo=try(maincMSGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
        cont=cont+1
        if(is.list(mo)){
            bicn=mo$ICL
            ICL[cont]=bicn}
        else{bicn=-Inf
            ICL[cont]=NA}
        if(bicn>bico){
            bico=bicn
            sg=G[b]
            model=mo
        }
    }
#    val=list(ICL=ICL,model=model)
    val=MixGHD(Index=ICL,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z, method="cMSGHD",data=as.data.frame(data),scale=scale)
    cat("The best model (ICL) for the range of  components used is  G = ", sg,".\nThe ICL for this model is ", bico,".",sep="")
    return(val)}
else if(modelSel=="AIC3"){
    bico=-Inf
    t=length(G)
    AIC3=matrix(NA,t,1)
    cont=0
    for(b in 1:t){
        mo=try(maincMSGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
        cont=cont+1
        if(is.list(mo)){
            bicn=mo$AIC3
            AIC3[cont]=bicn}
        else{bicn=-Inf
            AIC3[cont]=NA}
        if(bicn>bico){
            bico=bicn
            sg=G[b]
            model=mo
        }
    }
#    val=list(AIC3=AIC3,model=model)
    val=MixGHD(Index=AIC3,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z, method="cMSGHD",data=as.data.frame(data),scale=scale)
    cat("The best model (AIC3) for the range of  components used is  G = ", sg,".\nThe AIC3 for this model is ", bico,".",sep="")
    return(val)}
else {
    bico=-Inf
    t=length(G)
    AIC=matrix(NA,t,1)
    cont=0
    for(b in 1:t){
        mo=try(maincMSGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
        cont=cont+1
        if(is.list(mo)){
            bicn=mo$AIC
            AIC[cont]=bicn}
        else{bicn=-Inf
            AIC[cont]=NA}
        if(bicn>bico){
            bico=bicn
            sg=G[b]
            model=mo
        }
    }
    val=MixGHD(Index=AIC,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z, method="cMSGHD",data=as.data.frame(data),scale=scale)
    cat("The best model (AIC) for the range of  components used is  G = ", sg,".\nThe AIC for this model is ", bico,".",sep="")
    return(val)}






}

