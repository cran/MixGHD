      run_EMstep <- function(data=NULL, gpar=NULL, loglik=NULL, maxit=NULL, N=NULL, p=NULL, G=NULL, q=NULL, epsilon=NULL, label=NULL){
      counter =0
      mu = matrix(0, nrow=G, ncol=p)
      phi= matrix(0, nrow=G, ncol=p)
      alpha= matrix(0, nrow=G, ncol=p)
      cplr = matrix(0, nrow=G, ncol=p*2)
      gamr = matrix(0, nrow=G, ncol=p*p)
#      sigmar = matrix(0, nrow=G, ncol=p*p)
      cpl0 =matrix(0, nrow=G, ncol=2)
      wg = matrix(0, nrow=G, ncol=2)
      for(k in 1:G){
      par = gpar[[k]]
      mu[k,] = par$mu
      phi[k,] = par$phi
      alpha[k,] = par$alpha
      cplr[k,] = c(t(par$cpl))
      gamr[k,] = c(t(par$gam))
#      sigmar[k,] = c(t(par$sigma))
      cpl0[k,] = par$cpl0
      wg[k,] = par$wg
      }
      pi = gpar$pi
      v = 1
      mu1= t(mu)
      phi1=t(phi)
#      print("phi1")
#      print(phi1)
      alpha1=t(alpha)
#      print("alpha1")
#      print(alpha1)
      cpl1=t(cplr)
#      print("cplr")
#      print(cplr)
      gam1=t(gamr)
#      sigma1=t(sigmar)
      cpl01=t(cpl0)
      wg1=t(wg)
      if(is.null(label) ==T) label=rep(0,N)
      EM <- .C("EMstep", as.double(mu1), as.double(phi1), as.double(alpha1),
#            as.double(cpl1), as.double(gam1), as.double(sigma1), as.double(cpl01)
            as.double(cpl1), as.double(gam1), as.double(cpl01)
            , as.double(wg1), as.double(data), as.double(loglik), as.integer(maxit),
            as.integer(N), as.integer(p), as.integer(G), as.double(epsilon), as.integer(label),
            as.integer(v), as.double(pi), as.integer(counter), PACKAGE="MixGHD")
      loglik = EM[[9]]
      mu = matrix(EM[[1]], nrow=G,ncol=p, byrow=TRUE)
      phi = matrix(EM[[2]], nrow=G, ncol=p, byrow=TRUE)
      alpha = matrix(EM[[3]], nrow=G, ncol=p, byrow=TRUE)
      counter = EM[[18]]
      cpl = array(EM[[4]], dim=c(p,2,G))
      gam = array(EM[[5]], dim=c(p,p,G))
#      sigma = array(EM[[6]], dim=c(p,p,G))
      cpl0 = matrix(EM[[6]], nrow=G, ncol=2, byrow=TRUE)
      wg = matrix(EM[[7]], nrow=G, ncol=2, byrow=TRUE)
      gpar = list()
      for(k in 1:G){
         gpar[[k]] = list()
         gpar[[k]]$mu = mu[k,]
         gpar[[k]]$phi = phi[k,]
         gpar[[k]]$alpha = alpha[k,]
         gpar[[k]]$cpl = matrix(cpl[,,k],nrow=p, ncol=2,  byrow=TRUE)
         gpar[[k]]$gam = matrix(gam[,,k],nrow=p, ncol=p,  byrow=TRUE)
#         gpar[[k]]$sigma = matrix(sigma[,,k],nrow=p, ncol=p,  byrow=TRUE)
         gpar[[k]]$cpl0 = cpl0[k,]
         gpar[[k]]$wg = wg[k,]
         }
      gpar$pi = EM[[17]]
      val = list(loglik,gpar,counter)
      return(val)
     }
MainMCGHD=function(data=NULL, gpar0=NULL, G=2, max.iter=100, eps=1e-2,  label=NULL, method="km",nr=NULL){
    pcol=ncol(data)

#    if (is.null(gpar0)) gpar  = rmgpar(g=G,p=ncol(data),data=data, method=method,nr=nr)
#    else gpar = gpar0
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgparC(data=data, g=G, w=z,l=lc)

    }
    else{
    if (is.null(gpar0)) gpar  = rmgpar(g=G,p=ncol(data),data=data, method=method,nr=nr)
        else{ gpar = gpar0
            for(i in 1:G){
                gpar[[i]]$gam   = eigen( gpar0[[i]]$sigma)$vectors
                gpar[[i]]$phi   = eigen( gpar0[[i]]$sigma)$values}
    }

    }



    loglik = numeric(max.iter)
    for (i in 1:3) {
        gpar = EMgrstep(data=data, gpar=gpar, v=1, label = label,it=i)
        loglik[i] = llik(data, gpar)
    }
       N = nrow(data)
       p = ncol(data)
       maxit = max.iter
       temp <- run_EMstep(data, gpar, loglik, maxit, N, p, G, q, eps, label)
     i = temp[[3]]
#    if(i<max.iter){loglik=loglik[-(i+1:max.iter)]}
    if(i<max.iter){temp[[1]]=temp[[1]][-(i+1:max.iter)]}  
    BIC=2*temp[[1]][i]-log(nrow(data))*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    AIC=2*temp[[1]][i]-2*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    AIC3=2*temp[[1]][i]-3*((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2))
    z=weights(data=data, gpar= temp[[2]])
    ICL=BIC+2*sum(log(apply(z,1,max)))
    par=partrue(temp[[2]],G)
    val = list(loglik= temp[[1]][1:i], gpar=temp[[2]],par=par, z=z, map=MAP(data=data, gpar= temp[[2]], label=label),BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3 )
    return(val)

}





MCGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, eps=1e-2,  label=NULL, method="km",scale=TRUE,nr=10, modelSel="AIC" ) {
	data=as.matrix(data)
    if( scale==TRUE){
        data=scale(data)}
    	pcol=ncol(data)
        #if (nrow(data)<((G-1)+G*(4*pcol+2+pcol*(pcol-1)/2)))stop('G is too big, number of parameters > n')
	if (is.null(data)) stop('data is null')
	if (nrow(data) == 1) stop('nrow(data) is equal to 1')
	#if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
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
    mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
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
#########
#     val=list(BIC=BIC,model=model)
     val=list(index=BIC,model=model)
    cat("The best model (BIC) for the range of  components used is  G = ", sg,".\nThe BIC for this model is ", bico,".",sep="")
        return(val)}
   else if(modelSel=="ICL"){
        
        iclo=-Inf
        t=length(G)
        ICL=matrix(NA,t,1)
        cont=0
        for(b in 1:t){
            mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
            cont=cont+1
            if(is.list(mo)){
                icln=mo$ICL
                ICL[cont]=icln}
            else{icln=-Inf
                ICL[cont]=NA}
            if(icln>iclo){
               iclo=icln
                sg=G[b]
                model=mo
            }
        }
######        val=list(ICL=ICL,model=model)
        val=list(index=ICL,model=model)
        cat("The best model (ICL) for the range of  components used is  G = ", sg,".\nThe ICL for this model is ", iclo,".",sep="")
        return(val)}
   else if(modelSel=="AIC3"){
       
       iclo=-Inf
       t=length(G)
       AIC3=matrix(NA,t,1)
       cont=0
       for(b in 1:t){
           mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
           cont=cont+1
           if(is.list(mo)){
               icln=mo$AIC3
               AIC3[cont]=icln}
           else{icln=-Inf
               AIC3[cont]=NA}
           if(icln>iclo){
               iclo=icln
               sg=G[b]
               model=mo
           }
       }
#####       val=list(AIC3=AIC3,model=model)
       val=list(index=AIC3,model=model)
       cat("The best model (AIC3) for the range of  components used is  G = ", sg,".\nThe AIC3 for this model is ", iclo,".",sep="")
       return(val)}
   
   else {
       
       iclo=-Inf
       t=length(G)
       AIC=matrix(NA,t,1)
       cont=0
       for(b in 1:t){
           mo=try(MainMCGHD(data=data, gpar0=gpar0, G=G[b], max.iter, eps,  label, method,nr=nr),silent = TRUE)
           cont=cont+1
           if(is.list(mo)){
               icln=mo$AIC
               AIC[cont]=icln}
           else{icln=-Inf
               AIC[cont]=NA}
           if(icln>iclo){
               iclo=icln
               sg=G[b]
               model=mo
           }
       }
#####       val=list(AIC=AIC,model=model)
       val=list(index=AIC,model=model)
       cat("The best model (AIC) for the range of  components used is  G = ", sg,".\nThe AIC for this model is ", iclo,".",sep="")
       return(val)}
	}

