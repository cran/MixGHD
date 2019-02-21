      run_EMstepGH <- function(data= NULL, gpar = NULL, loglik=NULL, maxit=NULL, N= NULL, p=NULL, G=NULL, eps=NULL, label=NULL){
     counter =0
     mu = matrix(0,nrow=G,ncol=p)
     alpha = matrix(0,nrow=G,ncol=p)
     sigmar = matrix(0, nrow=G, ncol=p*p)
     cpl = matrix(0, nrow = G, ncol=2)
     for(k in 1:G){
       par = gpar[[k]]
       mu[k,] = par$mu
       alpha[k,] = par$alpha
       cpl[k,]= par$cpl
       sigmar[k,] = c(t(par$sigma))
    }
    pi = gpar$pi
    v=1
    mu1=t(mu)
    cpl1=t(cpl)
    alpha1=t(alpha)
    sigma1=t(sigmar)
    if(is.null(label) ==T) label=rep(0,N)
    EMGH <- .C("EMstepGH", as.double(mu1),as.double(alpha1),
                  as.double(sigma1),as.double(cpl1),as.integer(N), as.integer(p),
                  as.integer(G), as.double(loglik), as.integer(maxit), as.double(eps),as.integer(label), as.double(pi), as.integer(v), as.double(data), as.integer(counter), PACKAGE="MixGHD")
    loglik = EMGH[[8]]
    counter = EMGH[[15]]
    sigma = array(EMGH[[3]],dim=c(p,p,G))
    alpha = matrix(EMGH[[2]], nrow=G, ncol=p, byrow=TRUE)
    cpl = matrix(EMGH[[4]], nrow=G,ncol=2, byrow=TRUE)
    mu = matrix(EMGH[[1]],nrow =G, ncol=p, byrow=TRUE)
    gpar = list()
    for (k in 1:G){
        gpar[[k]] = list()
        gpar[[k]]$mu = mu[k,]
        gpar[[k]]$alpha = alpha[k,]
        gpar[[k]]$cpl = cpl[k,]
        gpar[[k]]$sigma = matrix(sigma[,,k], nrow=p, ncol=p, byrow=TRUE)
    }
    gpar$pi = EMGH[[12]]
    val = list(loglik,gpar,counter)
    return(val)

    }
mainMGHD<-function(data=NULL, gpar0, G, n, label  , eps, method  ,nr=NULL) {

    pcol=ncol(data)
    if(!is.null(label)){
        lc=apply(data[label==1,],2,mean)
       # if(min(label)==0&max(label)==G){
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        #}
        # else{
        #   print("G needs to be equal to max(label)")
        #   # for(i in 2:max(label)){
        #   #   lc=rbind(lc,apply(data[label==i,],2,mean))
        #   # }
        #   # for(i in (max(label)+1):G){
        #   #   lc=rbind(lc,apply(data[label==i,],2,mean))
        #   # }
        # }
       
        z = combinewk(weights=matrix(1/G,nrow=nrow(data),ncol=G), label=label)
        if (is.null(gpar0)) gpar  = rgparGH(data=data, g=G, w=z,l=lc)
        else gpar  = gpar0

    }
    else{
        if (is.null(gpar0)) gpar = try(igpar(data=data, g=G, method=method,nr=nr))
        else gpar  = gpar0}


    loglik = numeric(n)
    for (i in 1:3) {
        gpar = try(EMgrstepGH(data=data, gpar=gpar, v=1, label = label))        ###parameter estimation
        loglik[i] = llikGH(data, gpar)}
     maxit = n
     N = nrow(data)
     p = ncol(data)
 temp <- run_EMstepGH(data, gpar, loglik, maxit, N, p, G, eps, label)

    i = temp[[3]]
    if(i<n){temp[[1]]=temp[[1]][-(i+1:n)]}
    BIC=2*temp[[1]][i]-log(nrow(data))*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    z=weightsGH(data=data, gpar= temp[[2]])
    ICL=BIC+2*sum(log(apply(z,1,max)))
    AIC=2*temp[[1]][i]-2*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    AIC3=2*temp[[1]][i]-3*((G-1)+G*(2*pcol+2+pcol*(pcol-1)/2))
    val = list(loglik= temp[[1]], gpar=temp[[2]], z=z, map=MAPGH(data=data, gpar= temp[[2]], label=label),BIC=BIC,ICL=ICL,AIC=AIC,AIC3=AIC3 )
    return(val)

}


MGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL , eps=1e-2, method="kmeans" ,scale=TRUE ,nr=10, modelSel="AIC") {
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
	#if ( G < 1) stop('G is not a positive integer')
	if ( max.iter < 1) stop('max.iter is not a positive integer')
	
     if(modelSel=="BIC"){
    bico=-Inf
    t=length(G)
    BIC=matrix(NA,t,1)
    cont=0
	for(b in 1:t){
        mo=try(mainMGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
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
#    val=list(BIC=BIC,model=model)
    val=MixGHD(Index=BIC,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHD",data=as.data.frame(data),scale=scale)
    
    cat("The best model (BIC) for the range of  components used is  G = ", sg,".\nThe BIC for this model is ", bico,".",sep="")
         return(val)}
     
     
      else if(modelSel=="ICL"){
          bico=-Inf
          t=length(G)
          ICL=matrix(NA,t,1)
          cont=0
          for(b in 1:t){
              mo=try(mainMGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
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
#          val=list(ICL=ICL,model=model)
          val=MixGHD(Index=ICL,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHD",data=as.data.frame(data),scale=scale)
          
          cat("The best model (ICL) for the range of  components used is  G = ", sg,".\nThe ICL for this model is ", bico,".",sep="")
          return(val)}
      else if(modelSel=="AIC3"){
          bico=-Inf
          t=length(G)
          AIC3=matrix(NA,t,1)
          cont=0
          for(b in 1:t){
              mo=try(mainMGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
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
#          val=list(AIC3=AIC3,model=model)
          val=MixGHD(Index=AIC3,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHD",data=as.data.frame(data),scale=scale)
          
          cat("The best model (AIC3) for the range of  components used is  G = ", sg,".\nThe AIC3 for this model is ", bico,".",sep="")
          return(val)}
      else {
          bico=-Inf
          t=length(G)
          AIC=matrix(NA,t,1)
          cont=0
          for(b in 1:t){
              mo=try(mainMGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method,nr=nr),silent = TRUE)
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
          val=MixGHD(Index=AIC,AIC=model$AIC,AIC3=model$AIC3,BIC=model$BIC,ICL=model$ICL, map=model$map, gpar=model$gpar, loglik=model$loglik, z=model$z,method="MGHD",data=as.data.frame(data),scale=scale)
        
          #val=list(index=AIC,model=model)
          cat("The best model (AIC) for the range of  components used is  G = ", sg,".\nThe AIC for this model is ", bico,".",sep="")
          return(val)}
}



