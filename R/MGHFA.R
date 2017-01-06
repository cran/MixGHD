        run_EMstepFA <-function(data=NULL, gpar= NULL, loglik=NULL, maxit=NULL, N=NULL, p=NULL, G=NULL, q=NULL, epsilon=NULL, label=NULL){
        counter = 0
        mu = matrix(0,nrow = G, ncol=p)
        alpha = matrix(0, nrow = G, ncol=p)
        cpl = matrix(0, nrow = G, ncol=2)
        Lambdar = matrix(0, nrow = G, ncol=p*q)
        sigmar = matrix(0, nrow = G, ncol=p*p)
        errr = matrix(0, nrow=G, ncol =p*p)
        for(k in 1:G){
        par = gpar[[k]]
        mu[k,] = par$mu
        alpha[k,] = par$alpha
        cpl[k,]= par$cpl
        sigmar[k,] = c(t(par$sigma))
        Lambdar[k,] = c(t(par$Lambda))
        errr[k,]= c(t(par$err*(diag(p))))
        }
        pi = gpar$pi
        v=1
        mu1=t(mu)
        cpl1=t(cpl)
        alpha1=t(alpha)
        sigma1=t(sigmar)
        Lambda1=t(Lambdar)
        errr1=t(errr)
        if(is.null(label) ==T) label=rep(0,N)
        EMFA <- .C("EMstepFA",as.double(mu1),as.double(alpha1),
                  as.double(sigma1),as.double(cpl1),as.integer(N), as.integer(p),
                  as.integer(G), as.double(loglik), as.integer(maxit), as.double(epsilon), as.integer(label), as.double(pi),
                  as.double(Lambda1), as.double(errr1), as.integer(q), as.integer(v), as.double(data), as.integer(counter),PACKAGE="MixGHD")
        loglik = EMFA[[8]]
        counter = EMFA[[18]]
        sigma = array(EMFA[[3]],dim=c(p,p,G))
        alpha = matrix(EMFA[[2]], nrow=G, ncol=p, byrow=TRUE)
        cpl = matrix(EMFA[[4]], nrow=G,ncol=2, byrow=TRUE)
        mu = matrix(EMFA[[1]],nrow =G, ncol=p, byrow=TRUE)
        Lambda = array(EMFA[[13]], dim=c(p,q,G))
        err = array(EMFA[[14]], dim=c(p,p,G))
        gpar = list()
        for (k in 1:G){
           gpar[[k]] = list()
           gpar[[k]]$mu = mu[k,]
           gpar[[k]]$alpha = alpha[k,]
           gpar[[k]]$cpl = cpl[k,]
           gpar[[k]]$sigma = matrix(sigma[,,k], nrow=p, ncol=p, byrow=TRUE)
           gpar[[k]]$Lambda = matrix(Lambda[,,k], nrow=p, ncol=q, byrow=TRUE)
           gpar[[k]]$err = matrix(err[,,k], nrow=p, ncol=p, byrow=TRUE)
        }
        gpar$pi = EMFA[[12]]
        val = list(loglik,gpar,counter)
        return(val)

        }
mainMGHFA<-function(data=NULL, gpar0, G, n, label  , eps, method ,q,nr=nr ) {
    pcol=ncol(data)
    if(!is.null(label)&&min(label>0)){
        lc=apply(data[label==1,],2,mean)
        for(i in 2:G){
            lc=rbind(lc,apply(data[label==i,],2,mean))
        }
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=label)
        gpar  = rgpar(data=data, g=G, w=z,l=lc)

    }
    else{
        if (is.null(gpar0)) gpar = igparM(data=data, g=G,q=q,method=method,nr=nr)
        else gpar  = gpar0}

        loglik = numeric(n)
        for (i in 1:3) {
                gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label = label)     ###parameter estimation
                loglik[i] = llikFA(data, gpar) ##likelyhood
        }
        N = nrow(data)
        p = ncol(data)
        maxit = n
        temp <- run_EMstepFA(data, gpar, loglik, maxit, N, p, G, q, eps, label)



#    if(i<n){loglik=loglik[-(i+1:n)]}
     i=temp[[3]]
     if(i <n){temp[[1]]=temp[[1]][-(i+1:n)]}
        BIC=2*temp[[1]][i]-log(nrow(data))*((G-1)+G*(3*pcol+2+pcol*q-q*(q-1)/2))
        val = list(loglik= temp[[1]], gpar=temp[[2]], z=weightsFA(data=data, gpar= temp[[2]]), map=MAPFA(data=data, gpar= temp[[2]], label=label) , BIC=BIC)
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



