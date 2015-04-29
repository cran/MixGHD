.packageName<-'MixGHD'

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                        Coaleased                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################




###################################### Main ##########################################




EMgrstep <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL, vg=NULL) {
    
	G= length(gpar$pi);
	for (k in 1:G ) {
		if (is.null(w)) w = weights(data=data, gpar=gpar,v=v)
		if (!is.null(label)) w = combinewk(weights=w, label= label)
		vg = iweights(data=data, gpar=gpar[[k]],v=v)
		ng = apply(w,2,sum)
        
        #  wg=matrix((apply((w[,k]*vg),2,sum)/sum(w[,k]*vg)),nrow(data),2,1)
		gpar[[k]] = updatemaScplp(y=data, par=gpar[[k]], weights=w[,k],iweights=vg, alpha.known=NULL, v=v)#iweights=vg,
		gpar[[k]]$wg = apply((w[,k]*vg),2,sum)/ng[k]
		
	}
	gpar$pi = apply(w,2,mean)
	return(gpar)
}

#####Stopping criteria

getall <- function(loglik) {
	if (length(loglik) <3) stop("must have at least 3 likelihood values")
	n = length(loglik)
	lm1 = loglik[n]
	lm  = loglik[(n-1)]
	lm_1  = loglik[(n-2)]
	am = (lm1 - lm)/(lm - lm_1)
	lm1.Inf = lm + (lm1 - lm)/(1-am)
	val = lm1.Inf - lm
	if (is.nan(val)) val=0
	if (val < 0) val= 1
	return( val )
}

################ rotate the parameters


partrue<-function(gpar,G=2){
	gparT=gpar
	grpraT=list()
	for(i in 1:G){
		par=gpar[[i]]
		sort=sort.int(par$phi,decreasing=T,index.return=T)
		par$phi=sort$x
		par$gam=par$gam[,sort$ix]
		
		gam=gpar[[i]]$gam
		
		par$mu=par$mu%*%t(gam)
		par$alpha=par$alpha%*%t(gam)
        #par$Sigma=gam%*%diag(par$sigma)%*%t(gam)
		
		gparT[[i]]=par
	}
	return(gparT)
	
}




##################################density functions###################################

#####Mixture of GH and MSGH



ddmsghyp <- function(data,par,log=FALSE,invS=NULL){
	wg = par$wg
	
	val1 = ddghyp(y=data,par=par,log=log,invS=invS)
	val2 = dmsghyp(y=data,par=par,log=log)
	val = wg[1]*val1 + wg[2]*val2
	
	return(val)
}

ddmmsghyp <- function(data, gpar) {
	logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
	for (k in 1:length(gpar$pi) ) logz[,k] = ddmsghyp(data=data, par=gpar[[k]], log=F, invS=NULL)
	val = log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi))
	return(val)
}





##Multivariate GH rotated



ddghyp <- function(y, par, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
	x     = y %*%(par$gam)
    # x = y
    # mu    = c((par$gam)%*%par$mu)
	mu = par$mu
	if (length(par$sigma) ==1) sigma = matrix(par$sigma,1,1)
	else sigma = par$sigma
    # alpha = c((par$gam)%*%par$alpha)
	alpha= par$alpha
	
	if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
	d = length(mu)
	omega = exp(sum(log(par$cpl0[1])))
	lambda = par$cpl0[2];
	
	if (is.null(invS)) invS = ginv(sigma)
	
	pa = omega + alpha %*% invS %*% alpha
	mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
	kx = sqrt(mx*pa)
	xmu = sweep(x,2,mu)
	
	lvx = matrix(0, nrow=nrow(x), 4)
	lvx[,1] = (lambda - d/2)*log(kx)
	lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
	lvx[,3] = xmu %*% invS %*% alpha
	
	lv = numeric(6)
	lv[1] = -1/2*log(det( sigma )) -d/2*(log(2)+log(pi))
	lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
	lv[3] = -lambda/2*( log(1) )
	lv[4] = lambda*log(omega) *0
	lv[5] = (d/2 - lambda)*log( pa )
	
	val = apply(lvx,1,sum) + sum(lv)
	if (!log) val = exp( val )
	
	return(val)
}

### MS GH


dmsghyp <- function(y, par, log=FALSE) {
    # x is a n x p matrix
	x     = y %*% (par$gam)
	mu    = par$mu; phi = par$phi;
	alpha = par$alpha;
	d = length(mu); chi = par$cpl[,1]; psi = par$cpl[,1];
	lambda = par$cpl[,2];
	
	xmu = sweep(x,2,mu,"-")
	
	pa = psi + alpha^2/phi   # p numeric
	mx = sweep(sweep(xmu^2,2,1/phi, FUN="*"), 2,chi, "+") # n x p matrix
	kx = sqrt(sweep(mx, 2, pa, "*")) # nxp matrix
	
	lx1 = sweep( sweep(log(mx),2,log(pa),"-"), 2, (lambda - 1/2)/2, "*")
	lx2 = t(apply(kx, 1, function(z,lam=NULL) { log(besselK( z, nu=lambda-1/2, expon.scaled =TRUE)) - z }, lam=lambda ))
	lx3 = sweep(xmu, 2, alpha/phi, FUN="*")
	
	lv = numeric(3)
	lv1 = -1/2*log( phi ) -1/2*(log(2)+log(pi)) # d=1
	lv2 = - (log(besselK( sqrt(chi*psi), nu=lambda, expon.scaled =TRUE)) - sqrt(chi*psi) )
	lv3 = lambda/2*(log(psi)-log(chi) )
	
	if(ncol(y)==1){lx2=t(lx2)}
    lx=lx1 +lx2 + lx3
    #val = apply(lx1 + lx2 + lx3, 1, sum) + sum(lv1 + lv2 + lv3)
    val = apply(lx, 1, sum) + sum(lv1 + lv2 + lv3)
    
	if (!log) val = exp( val )
	
	return(val)
}

##################################### Log likelihood ##################################



llik <- function(data, gpar) {
	logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
	for (k in 1:length(gpar$pi) ) logz[,k] = ddmsghyp(data=data, par=gpar[[k]], log=TRUE, invS=NULL)
	val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
	return(val)
}







####################################labeling##########################################

MAP <- function(data, gpar, label=NULL) {
	w = weights(data=data, gpar=gpar, v=1)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
	z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
	z = as.numeric(z)
	return( z)
}


######################################### E step #######################################

##############zig  model weights

combinewk <- function(weights=NULL, label=NULL)	{
    # known is a numeric with
    # 0 if unknown group membership
    # 1,2,3,.. for label of known group
	if (is.null(label)) stop('label is null')
	kw     = label !=0
	for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
	return(weights)
}



weights <- function(data=NULL, gpar=NULL, v=1) {
	G = length(gpar$pi)
	if (G > 1) {
		zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
		for (k in 1:G ) zlog[,k] =  ddmsghyp(data=data, par=gpar[[k]], log=TRUE)
		w = t(apply(zlog, 1, function(z,wt,v) {
            x=exp(v*(z + log(wt)) );
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
	} else w = matrix(1,nrow=nrow(data), ncol=G)
	return(w)
}


##############uig  mixtures weights

iweights <- function(data=NULL, gpar=NULL, v=1) {
	#G = 2
	#if (G > 1) {
    zlog = matrix(0, nrow=nrow(data), ncol=2)
    zlog[,1] =  ddghyp(y=data, par=gpar, log=TRUE, invS=NULL)
    zlog[,2] =  dmsghyp(y=data, par=gpar, log=TRUE)
    w = t(apply(zlog, 1, function(z,wt,v) {
        x=exp(v*(z + log(wt)) );
        
        if (sum(x)  == 0) x= rep(1,length(x))
        x =  x/sum(x)
        return( x )
    }, wt=gpar$wg,v=v ))
    #} else w = matrix(1,nrow=nrow(data), ncol=G)
	return(w)
}


############Expected values GIG


###Bessel function
logbesselKv <- function(x, y) { log(besselK(x=y, nu=x, expon.scaled=TRUE)) - log(y)}
logbesselKvFA <- function(x, y) {
    val = log(besselK(x=y, nu=x, expon.scaled=FALSE)) - y
    sun = is.infinite(val)
    val[sun] = besselK.nuAsym(x=y[sun], nu=abs(x[sun]), k.max=4, log=TRUE, expon.scaled=FALSE) - y[sun]
    return(val)
}
besselKv    <- function(x, y) { besselK(x=y, nu=x)}


######GIG univariate
gig <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3
	sab =  sqrt(a*b)
	kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
	kv  = besselK( sab, nu=v, expon.scaled =TRUE)
	kv12 = kv1/kv
	
	sb.a = sqrt(b/a)
	w    = kv12*sb.a
	invw = kv12*1/sb.a - 2*v/b
	logw = log(sb.a) + grad( logbesselKvFA, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
	
	
	val =cbind(w,invw,logw)
	return(val)
}

gig2 <- function(x=NULL, par=NULL, invS=NULL) {
    # returns the same as gig
	d = length(par$mu)
	
	if (is.null(invS)) invS = ginv(par$sigma)
	alpha = par$alpha
	omega = exp(  log(par$cpl[1])  )
	
	a1 = omega + as.numeric( alpha %*% invS %*% alpha )
	b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
	v1 = par$cpl0[2]-d/2
	
	val = gig(b=b1,a=a1, v=v1)
	return(val)
}


#####GIG MS
gigp <- function(a=NULL,B =NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3
	SaB =  sqrt(sweep(B, 2, a, "*" ))
	
	Kv12 = t(apply(SaB, 1, function(x=NULL, v=NULL) {
        kv1 = besselK( x, nu=v+1, expon.scaled =TRUE)
        kv0 = besselK( x, nu=v, expon.scaled =TRUE)
        kv12 = kv1/kv0
        return(kv12)
    }, v=v ) )
	
	SB.a =  sqrt(sweep(B, 2, 1/a, "*" ))
    if(nrow(Kv12)==1){Kv12=t(Kv12)}
	W    = Kv12*SB.a
	invW = Kv12/SB.a - 2*sweep(1/B, 2, v, "*" )
    #print(cbind(rep(v,each=nrow(SaB)),as.numeric(SaB)) )
	
	logW = log(SB.a) + matrix(grad( logbesselKvFA, x=rep(v,each=nrow(SaB)), y=as.numeric(SaB), method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)), nrow=nrow(SaB), ncol=ncol(SaB) )
	
	val = list(W=W,invW=invW,logW=logW)
	return(val)
}


gig2p <- function(sr2=NULL, par=NULL) {
    # returns the same as gig
	
	omega = exp( (log(par$cpl[,1])) )
	a1 = omega + par$alpha^2/par$phi   # vector
	B1 = sweep(sr2, 2, omega, "+") # matrix
	v1 = par$cpl[,2]-1/2 # vector
	
	val = gigp(a=a1, B=B1, v=v1)
	return(val)
}



weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )

########################### M step Parameter estimation ##################################

################Initialization
rmgpar <- function(g=NULL, p=NULL,data=NULL, method="kmeans",n=10) {
	val = list()
	#set.seed(142857)
    if(method=="modelBased"){
        l=gpcm(data,  G=g, mnames=c("VVV"))$map
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z)
        
        for (j in 1:n)  gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)
        # l=BAR(data,l)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z)
        
        for (j in 1:n)  gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z)
        return(gpar)
        #l=BAR(data,l)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        l=round(runif(nrow(data))*(g-1)+1)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z)
        #l=BAR(data,l)
    }
    else{ l=kmeans(data,g)
        l=l$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparC(data=data, g=g, w=z)
        
        for (j in 1:n)  gpar = EMgrstep(data=data, gpar=gpar, v=1, label= l,  w=z)}
    return(gpar)
    
}

rgparC<- function(data,g=NULL, w=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
	val = list()
    # l=rmgparMS(g, p,data,method=method)
	for (k in 1:g) val[[k]] = rmparC(data=data,wt=w[,k],k)
	val$pi = rep(1/g,g)
	return(val)
}

rmparC <- function(data,wt,k) {
    par = list()
    p=ncol(data)
	par$mu =rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )#l[k,] #rnorm(p,0,.01);
	#par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
	par$alpha =  rnorm(p, 0, sqrt(1/nrow(data)) )#rnorm(p,0,.01);
    
    sigma = ( ( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
    par$sigma=sigma
	par$cpl0 = c(1,-1/2)
    par$cpl = cbind( rep(1,p), rep(-1/2,p))
	wg <- runif(1)
	par$wg = c(wg,1-wg)
	return(par)
}
########################## MAIN

updatemaScplp <- function(y=NULL, par=NULL, weights=NULL, iweights=NULL, alpha.known=NULL, v=1) {
	p=ncol(y)
	n=nrow(y)
	if (is.null(weights)) weights=rep(1,n)
	if (is.null(iweights)) iweights=matrix(c(0.5),nrow=n,ncol=2)
	
	x = y %*% (par$gam)
	sr2 = sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
	W <- matrix(0,n,p)
	invW <- matrix(0,n,p)
	logW <- matrix(0,n,p)
	abcM <- list(W=W,invW=invW,logW=logW)
    #	A <- 0
    if(sum(iweights[,2]) != 0){
        abcM = gig2p(sr2=sr2, par=par)
        #		A = apply(abcM$W,   2,weighted.sum, wt=weights*iweights[,2])
        #		B = apply(abcM$invW,2,weighted.sum, wt=weights*iweights[,2])/sum(weights*iweights[,2])
        #		C = apply(abcM$logW,2,weighted.sum, wt=weights*iweights[,2])/sum(weights*iweights[,2])
    }
	
	abc = gig2(x=x, par=par)
	
	if (runif(1) > 1/2) new.gam = updategam2(gam0=par$gam, y=y, phi=par$phi,alpha=par$alpha, mu=par$mu, wt=weights,ut=iweights, invW= abcM$invW,cpl0=abc)
	else new.gam = updategam1(gam0=par$gam, y=y,  phi=par$phi,alpha=par$alpha, mu=par$mu, wt=weights,ut=iweights, invW=abcM$invW,cpl0=abc)
	
    
    
    x = y %*% (new.gam)
	
	sumw = sum(weights)
  	
    #	AU = sum((abc[,1]*weights*iweights[,1]))
    #	BU = sum((abc[,2]*weights*iweights[,1]))
    #	CU = sum((abc[,3]*weights*iweights[,1]))
	#Z.x <- apply(x,2,weighted.sum,wt=weights*iweights[,1])
	
    #############MU
    S1=abcM$W*iweights[,2]+abc[,1]*iweights[,1]
    S2=abcM$invW*iweights[,2]+abc[,2]*iweights[,1]
    S=S1*S2
	sS1=apply(S1,2,weighted.sum, wt=weights)#(apply((abcM$W*weights*iweights[,2]),2,sum)+sum(abc[,1]*weights*iweights[,1]))
	sS2=apply(S2,2,sum)#(apply((abcM$invW*weights*iweights[,2]),2,sum)+sum(abc[,2]*weights*iweights[,1]))#/sumw
    sumS=apply(S,2,weighted.sum, wt=weights^2)#/sumw
    #molt=(sS1*S2-sumw)
    #moltb=(sS1*sS2-sumw^2)
    molt=(S1*S2-1)#S1*S2
    moltb=(sumS-sumw^2)
    nummu=apply(weights^2*x*molt,2,sum)
    #denmu=apply(weights*moltb,2,sum)
    mu.new=nummu/moltb
    mu.new=as.vector(mu.new)
	
    #	S2=apply((abcM$invW*weights*iweights[,2]),2,sum)+rep(sum(abc[,2]*weights*iweights[,1]),1,p)
    #mu.num <- t(apply(x, 2, weighted.sum, w=abcM$invW*weights*iweights[,2]) +apply(x, 2, weighted.sum, w=abc[,2]*weights*iweights[,1])-#sum(weights)*ginv(S1)*apply(x, 2, weighted.sum, w=weights))
    #	mu.den <- sum(abcM$invW*weights*iweights[,2]) + sum(abc[,2]*weights*iweights[,1])-sumw^2*ginv(S1)
    #	mu.new = t(mu.num)/mu.den
    #	mu.new=as.vector(mu.new)
    # #######ALPHA
    #sumS2=apply(S2,2,weighted.sum, wt=weights)
    
    #   molt2=sumS2/sumw-S2
    #   numal=apply(weights*x*molt2,2,sum)
    #alpha.new=numal/denmu
    
	
    Numal=apply((sweep(x, 2, mu.new,"-")), 2,weighted.sum, wt=weights)
	alpha.new=Numal/sS1
    
    
    #########PHI
	
    #	Num11 = 0
    #	if((sum(iweights[,1]) != 0) & sum(abc[,2]*weights*iweights[,1])!=0) Num11 = cov.wt(x=x, wt=(abc[,2]*weights*iweights[,1]), center=mu.new, method="ML")$cov*BU
    #	r = as.numeric(Z.x - mu.new)
	#m = as.numeric(alpha.new)
    #	Num12 = outer(alpha.new,r)
    #	Num13 = outer(r, alpha.new)
	#Num14 = outer(alpha.new, alpha.new)*AU
	
	#Ax = apply(sweep(x, 2, mu.new,"-")^2*abcM$invW, 2, weighted.sum, wt=weights*iweights[,2])
	#ax = apply(x, 2, weighted.sum, wt=weights*iweights[,2])
	#phi.new1 = (Num11 - (Num12 + Num13) + Num14 + Ax - (2*(ax - mu.new)*alpha.new + alpha.new^2*A))/(sum(weights*iweights[,1]) + sum(weights*iweights[,2]))
	
    # phi.new = diag(phi.new1)
	
    #phi.new=par$phi
	
	
    #	InvM=matrix(abc[,2],nrow(x),p)
    # update for sigma
    #	Matu = x-matrix(mu.new,n,p,1)-matrix(abc[,1],n,p)*matrix(alpha.new,n,p,1)
    #	MatM = x-matrix(mu.new,n,p,1)-abcM$W*matrix(alpha.new,n,p,1)
	
	#Num1sig=t(weights*iweights[,1]*abc[,2]*Matu)%*%(Matu)
	
	#Num2sig=t(weights*iweights[,2]*abcM$invW*MatM)%*%(MatM)
    
    r=x-matrix(mu.new,n,p,1)
    Num1sig=t(matrix(weights*iweights[,1]*abc[,2],n,p)*r)%*%r+(sum(weights*iweights[,1]*abc[,1])*alpha.new)%*%t(alpha.new)-(apply(matrix(weights*iweights[,1],n,p)*r,2,sum))%*%t(alpha.new)-(alpha.new)%*%t(apply(matrix(weights*iweights[,1],n,p)*r,2,sum))
    
    
    Num2sig=t(matrix(weights*iweights[,2],n,p)*abcM$invW*r)%*%r+(apply(weights*iweights[,2]*abcM$W,2,sum)*alpha.new)%*%t(alpha.new)-(apply(matrix(weights*iweights[,2],n,p)*r,2,sum))%*%t(alpha.new)-(alpha.new)%*%t(apply(matrix(weights*iweights[,2],n,p)*r,2,sum))
	
    phi.new=diag(Num1sig+Num2sig)/sumw
    
    #Ax = apply(sweep(x, 2, mu.new, "-")^2*abcM$invW, 2, weighted.sum, wt=weights*uwe)/sumw
    #	ax = apply(x, 2, weighted.sum, wt=weights*uwe)/sumw
    #	AxU = apply(sweep(x, 2, mu.new, "-")^2*InvM, 2, weighted.sum, wt=weights*uweU)/sumw
    #	axU = apply(x, 2, weighted.sum, wt=weights*uweU)/sumw
    #	sigma.new = Ax - 2*(ax - mu.new)*alpha.new + alpha.new^2*A+AxU - 2*(axU - mu.new)*alpha.new + alpha.new^2*AU
    #sigma.new=par$sigma
    
    ################CPLMs
    if(sum(iweights[,2]) != 0){
        abcM = gig2p(sr2=sr2, par=par)
        A = apply(abcM$W,   2,weighted.sum, wt=weights)/sum(weights)
        B = apply(abcM$invW,2,weighted.sum, wt=weights)/sum(weights)
        C = apply(abcM$logW,2,weighted.sum, wt=weights)/sum(weights)
        
		
		omega=  (par$cpl[,1])
		test  = cbind(omega, lambda=par$cpl[,2], A, B, C)
        # print(test)
		cpl.new = t(apply(test, 1, function(z) {
            temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
            return( c(  temp[1], temp[2]) )
        }))
  	} else {
  		cpl.new=par$cpl
	}
    ##############CPL uni
	
	
    
	par.ol = par$cpl0
	ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    
	a = updateol(ol=par.ol, ABC=ABC, n=2)
	cpl.newU =  c( a[1],  a[2])
	
    # cpl.newU=par$cpl0
    
    #	cpl.new=par$cpl
    # phi.new=par$phi
    # alpha.new=par$alpha
    #    mu.new=par$mu
    #   new.gam=par$gam
    par$sigma=new.gam%*%(diag(p)*phi.new)%*%t(new.gam)
	
	new.par = list(mu=mu.new, phi=phi.new, alpha=alpha.new, cpl=cpl.new, gam=new.gam, sigma=par$sigma, cpl0=cpl.newU )
	return(new.par)
}



######################## MS cpl






updateol <- function(ol=NULL, ABC=NULL, n=1) {
	ol0 = ol
	for (i in 1:n) {
		if (ABC[3] == 0) {
			ol[2] = 0
		} else {
            #print( c(1, ol[2], ol[1] ) )
            #if(ol[1]<=0){ol[1]=eps}
			bv = grad( logbesselKvFA, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
			ol[2] = ABC[3]*(ol[2]/bv)
		}
		
		Rp = Rlam(ol[1],lam=+ol[2])
		Rn = Rlam(ol[1],lam=-ol[2])
		f1 = Rp + Rn - (ABC[1]+ABC[2])
		f2 = ( Rp^2 - (2*ol[2]+1)/ol[1] *Rp -1 ) + ( Rn^2 - (2*(-1*ol[2])+1)/ol[1] *Rn -1 )
		ol=round(ol,4)
        
        if ( ol[1] - f1/f2 > 0 ) ol[1] = ol[1] - f1/f2
        #ol[1] =  ol[1] - f1/f2
        
    }
    return(ol)
}


Rlam <- function(x, lam=NULL) {
	v1 = besselK(x, nu= lam+1, expon.scaled=FALSE)
	v0 = besselK(x, nu= lam, expon.scaled=FALSE)
	val = v1/v0
    
    set= is.infinite(v1)|is.infinite(v0)|v0==0|v1==0
    if(any(set)){
        lv1=besselK.nuAsym(x=x[set],nu=abs(lam[set]+1),k.max=4, expon.scaled=FALSE,log=TRUE)
        lv0=besselK.nuAsym(x=x[set],nu=abs(lam[set]),k.max=4, expon.scaled=FALSE,log=TRUE)
        val[set]=exp(lv1-lv0)
    }
	return(val)
}





################################## gamma





updategam1 <- function(gam0=NULL, y=NULL, phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL,cpl0=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/phi, "*" )
    Bs2 = matrix(cpl0[,2],nrow(x),ncol(x))/phi
    
    # u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt*ut[,2], "*")
    wty2 = sweep(y,1,wt*ut[,1],"*")
    
    F0t = t( x * Bs ) %*% wty + t(x*Bs2) %*% wty2
    e2 = apply(Bs,1,max)*wt*ut[,2]
    e1 = max(Bs2)*wt*ut[,1]
    if(sum(e2)==0){    F1t =  ((cov.wt(y, center=rep(0,ncol(y)), wt=e1, method="ML" )$cov*sum(e1)) %*% (gam0))
        
    } else {
        if(sum(e1)==0){F1t = ((cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0))
            
        } else {
            F1t = ((cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0)) + ((cov.wt(y, center=rep(0,ncol(y)), wt=e1, method="ML" )$cov*sum(e1)) %*% (gam0))
            
        }
    }
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    v1 = (Bs2*mu + alpha/phi)
    A0 = t(v)  %*% wty
    A1 = t(v1) %*% wty2
    
    
    F2 = F0t - t(F1t) - A0 - A1
    z  = svd( F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    #gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    if (obj.gam(gam0=gam,y=y, phi=phi, alpha=alpha, mu=mu,  wt=wt, ut=ut, invW=invW,abc=cpl0) <  obj.gam(gam0=gam0,y=y, phi=phi, alpha=alpha, mu=mu, wt=wt, ut=ut, invW=invW,abc=cpl0) ) {
        #		print('not good enough')
        gam = gam0
    }
    
    return(gam)
}



updategam2 <- function(gam0=NULL, y=NULL,  phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL,cpl0=NULL) {
    
    x = y %*% (gam0)
    sumw = sum(wt)
    Bs = sweep(invW, 2, 1/phi, "*" )
    Bs2 =  matrix(cpl0[,2],nrow(x),ncol(x))/phi
    # u  = sweep(y * invW, 1, wt, "*")
    wty = sweep(y, 1, wt*ut[,2], "*")
    wty2 = sweep(y,1,wt*ut[,1],"*")
    F0t = t( x * Bs ) %*% wty + t(x*Bs2) %*% wty2
    
    e1 = apply(y^2, 1, sum)*wt*ut[,1]
    e2 = apply(y^2, 1, sum)*wt*ut[,2]
    
    if(ncol(y)==1){F1t = (apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0) + weighted.sum(Bs2,wt=e1)*t(gam0)}
    else{
        F1t = diag(apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0) + weighted.sum(Bs2,wt=e1)*t(gam0)}
    
    
    #v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+") + (Bs2*mu + alpha/phi)
    #A0 = t(v)  %*% wty
    
    v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
    v1 = (Bs2*mu + alpha/phi)
    A0 = t(v)  %*% wty
    A1 = t(v1) %*% wty2
    
    
    F2 = F0t - F1t - A0 -A1
    z  = svd( F2/sumw)
    gam = ( (z$v) %*% t(z$u) )
    gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
    
    if (obj.gam(gam0=gam,y=y, phi=phi, alpha=alpha, mu=mu,  wt=wt, ut=ut, invW=invW,abc=cpl0) <  obj.gam(gam0=gam0,y=y, phi=phi, alpha=alpha, mu=mu, wt=wt, ut=ut, invW=invW,abc=cpl0) ) {
        #		print('not good enough')
        gam = gam0
    }
    
    return(gam)
}
tr <- function(A=NULL) { sum(diag(A)) }



obj.gam <- function(gam0=NULL, y=NULL, phi=NULL, alpha=NULL, mu=NULL, wt=NULL, ut=NULL, invW=NULL,abc=NULL) {
	
	invUmat=matrix(abc[,2],nrow(y),ncol(y))
    ###MS part
	x = y %*% (gam0)
	sumw = sum(wt)
	Bs  = sweep(invW, 2, 1/phi, "*" )
	wtx = sweep(x, 1, (wt*ut[,2]), "*")
	F0t = t( x * Bs ) %*% wtx
	
	v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
	A0 = t(v)  %*% wtx
	
	valMS = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
	
	
    ###Multivariate part
	x = y %*% (gam0)
	sumw = sum(wt)
	Bs  = sweep(invUmat, 2, 1/phi, "*" )
	wtx = sweep(x, 1, (wt*ut[,1]), "*")
	F0t = t( x * Bs ) %*% wtx
	
	v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/phi, "+")
	A0 = t(v)  %*% wtx
	
	valMUL = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
	
	val=valMS+valMUL
	
	return(val)
}





#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                       Factor analyzer                                                     ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



####################### Univ CPL



RRlamz <- function(x,lam=NULL,z=1) {
    val =RlamU(x, lam=-lam)+ RlamU(x, lam= lam)
    zval = val - z
    return(zval)
}

RlamU <- function(x, lam=NULL) {
    
    v1 = besselK(x, nu= lam+1, expon.scaled=TRUE)
    v0 = besselK(x, nu= lam, expon.scaled=TRUE)
    val = v1/v0
    return(val)
}

updateolU <- function(ol=NULL, ABC=NULL, n=1) {
    
    for (i in 1:n) {
        if (ABC[3] == 0) {
            ol[2] = 0
        } else {
            # if(ol[1]<=0){ol[1]=eps}
            bv = grad( logbesselKv, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
            ol[2] = ABC[3]*(ol[2]/bv)
        }
        
        
        lam0 = ol[2]
        omg0 = ol[1]
        Rp = RlamU(omg0,lam=+lam0)
        Rn = RlamU(omg0,lam=-lam0)
        f1 = Rp + Rn - (ABC[1]+ABC[2])
        f2 = ( Rp^2 - (2*lam0+1)/omg0 *Rp -1 ) + ( Rn^2 - (2*(-1*lam0)+1)/omg0 *Rn -1 )
        # note, it is suppose to be f1/2 and f2/2
        if ( ol[1] - f1/f2 > 0 ) ol[1] = ol[1] - f1/f2
    }
    
    return(ol)
}





weightsFA <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G ) zlog[,k] =  ddghypFA(x=data, par=gpar[[k]], log=TRUE)
        w = t(apply(zlog, 1, function(z,wt,v) {
            x=exp(v*(z + log(wt)) );
            sun=is.infinite(x)
            x[sun]=1
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}



###Function1




igparM <- function(data=NULL, g=NULL,q=2,method="kmeans") {
    ##initialization
    gpar = igpar3M(data=data, G=g, n=10,q=q,method=method)
    return(gpar)
}


####Function for igparM
igpar3M <- function(data=NULL, G=NULL, n=10,label=NULL,q=2,method="kmeans") {
    #set.seed(142857)
    if(method=="modelBased"){
        l=gpcm(data,  G=G, mnames=c("VVV"))$map}
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=G))
    }
    else{#set.seed(142857)
        l=kmeans(data,G)$cluster}
    
    lw = l#kmeans(data, centers=G, iter.max=10)$cluster
    z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
    gpar  = rgpar(data=data, g=G, w=z,q=q)
    
    for (j in 1:n) try({ gpar = EMgrstepFA(data=data, gpar=gpar, v=1, label= label,  w=z)}, TRUE)
    return(gpar)
}


rgpar <- function(data, g=2, w=NULL,q=2) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    for (k in 1:g) val[[k]] = iparM(data=data, wt=w[,k],q=q)
    val$pi = rep(1/g,g)
    return(val)
}


iparM <- function(data, wt=NULL,q=2) {
    if (is.null(wt)) wt = rep(1,nrow(data))
    p = ncol(data)
    val = list()
    val$mu    = rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )
    val$alpha = rnorm(p, 0, sqrt(1/nrow(data)) )
    e=eigen(cov.wt(data, wt = wt, method="ML")$cov)
    dia=diag(e$values[1:q])
    if(q==1){val$Lambda=as.matrix(e$vectors[,1],p,q)
        dia=e$vectors[,1:q]
        val$sigma =val$Lambda%*%t(val$Lambda)
    }
    else{
        val$Lambda=e$vectors[,1:q]
        val$sigma =val$Lambda%*%dia%*%t(val$Lambda)######Cris###################################################################################
    }
    val$err=diag(cov.wt(data, wt = wt, method="ML")$cov-val$sigma)*diag(p)########Cris########################################################################################################################
    #  val$sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    if (any(eigen(val$sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    val$cpl   = c(1,-1/2)
    return(val)
    
}
#####end function for igpar

EMgrstepFA <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
    ##parameter estimtion
    if (is.null(w)) w = weightsFA(data=data, gpar=gpar,v=v)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);##number clusters
    d= length(gpar[[1]]$mu);##number variables
    #Sk = array(0, c(d,d,G) )
    for (k in 1:G ) {
        
        gpar[[k]] = updatemaScplM1(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        w = weightsFA(data=data, gpar=gpar,v=v)
        gpar[[k]] = updatemaScplM2(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        #Sk[,,k]   = gpar[[k]]$sigma
        
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}








##################################################NEW###############################################################################

##################################################NEW###############################################################################






updatemaScplM1 <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2FA(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        A = ABC[1]
        B = ABC[2]
        u = (B - abc[,2])*weights
        t = (A*abc[,2]-1)*weights
        T = sum(t)
        
        mu.new    = apply(x, 2, weighted.sum, w=t)/T
        alpha.new = apply(x, 2, weighted.sum, w=u)/T
    } else {
        alpha.new = alpha.known
        mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
    }
    alpha.new = alpha.new*v
    
    
    
    
    par.old = c(  log(par$cpl[1]) , par$cpl[2])
    #a = optim( par.old, loglikgig4, ABC=ABC)
    par.ol = c(exp(par.old[1]), par.old[2])
    a = updateolU(ol=par.ol, ABC=ABC, n=2)
    #if ( loglikgig4(par.old, ABC) < loglikgig4(par.old, ABC) ) cpl.new = par$cpl
    #  else cpl.new =  c( rep(exp(a$par[1]), 2), a$par[2])
    #if(a[1]<=0){a[1]=eps}
    cpl.new =  c( a[1], a[2])
    
    
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=par$sigma, cpl=cpl.new, Lambda=par$Lambda, err=par$err )
    
    return(new.par)
}












updatemaScplM2 <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2FA(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    mu.new   = par$mu
    
    alpha.new = par$alpha
    
    
    cpl.new =  par$cpl
    
    
    
    
    #####second estep
    d=ncol(x)
    q=ncol(par$Lambda)
    var=par$Lambda
    fi=par$err
    fi1=ginv(fi)
    fi1=diag(d)*as.vector(fi1)
    dia=diag(q)
    if( length(q)==0){dia=q}
    inv=fi1-fi1%*%var%*%ginv(dia+t(var)%*%fi1%*%var)%*%t(var)%*%fi1
    beta=t(var)%*%inv
    term1=sweep(x,2, FUN ="-", STATS=mu.new)-matrix(alpha.new,nrow(x),d,1)*matrix(abc[,1],nrow(x),d)
    uhat=beta%*%t(term1)
    term1=sweep(x,2, FUN ="-", STATS=mu.new)*matrix(abc[,2],nrow(x),d)-matrix(alpha.new,nrow(x),d,1)
    
    uhatb=beta%*%t(sweep(x,2, FUN ="-", STATS=mu.new)*matrix(abc[,2],nrow(x),d)-matrix(alpha.new,nrow(x),d,1))
    
    
    A = cov.wt(x, wt=abs(abc[,2]*weights), center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data
    r = apply(x, 2, weighted.sum, wt=weights)/sumw
    R = A - (outer( r - mu.new, alpha.new) + outer(alpha.new, r - mu.new)) + outer(alpha.new,alpha.new)*ABC[1]
    
    
    euu=sum(abc[,2]*weights)*diag(q)-sum(abc[,2]*weights)*beta%*%var+beta%*%R%*%t(beta)*sumw
    
    #euu=diag(q)-beta%*%var+beta%*%t(term1)%*%term1%*%t(beta)
    va.new=(t(sweep(x,2, FUN ="-", STATS=mu.new)*(weights))%*%t(uhatb)-t(matrix(alpha.new,nrow(x),d,1)*(weights))%*%t(uhat))%*%ginv(euu)#/sum(abc[,2]*weights)
    
    # va.new=(t(sweep(x,2, FUN ="-", STATS=mu.new)*weights)%*%t(uhatb)-(alpha.new)%*%t(apply(weights*uhat,1,"sum")))%*%ginv(euu)
    
    
	#euu=diag(q)-beta%*%var+beta%*%t(term1*abc[,2]*weights)%*%term1%*%t(beta)
    fip1=t(sweep(x,2, FUN ="-", STATS=mu.new)*weights)%*%t(uhatb)%*%t(var)
    fip2=(alpha.new)%*%t(apply(weights*uhat,1,"sum"))%*%t(var)+(var%*%euu%*%t(var))
    
    fi.new=diag(R+(2/sumw)*(-fip1+fip2))
    # fi.new=diag(R+(t(-2*t(sweep(x,2, FUN ="-", STATS=mu.new)*weights))%*%t(uhatb)%*%t(var)+2*(alpha.new)%*%t(apply(weights*uhat,1,"sum"))%*%t(var)+(var%*%euu%*%t(var)))/sumw)
    # euu=diag(q)-beta%*%var+beta%*%t(term1*abc[,2]*weights)%*%term1%*%t(beta)
    #fi.new=diag(R-(t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)-sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    #fi.new=diag((t(abc[,2]*term1*weights)%*%(term1)-t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)+sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    # fi.new=1/fi.new
	#fi.new=par$err
	
    #va.new=R%*%t(inv)%*%var
    # inv=fi1-fi1%*%va.new%*%ginv(dia+t(va.new)%*%fi1%*%va.new)%*%t(va.new)%*%fi1
    #fi.new=(R%*%t(inv))*diag(d)%*%fi
    
    var=va.new%*%t(va.new)+diag(d)*fi.new
    
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=var, cpl=cpl.new, Lambda=va.new, err=fi.new )
    
    return(new.par)
}







updatemaScplMbu <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL,pi=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2FA(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        A = ABC[1]
        B = ABC[2]
        u = (B - abc[,2])*weights
        t = (A*abc[,2]-1)*weights
        T = sum(t)
        
        mu.new    = apply(x, 2, weighted.sum, w=t)/T
        alpha.new = apply(x, 2, weighted.sum, w=u)/T
    } else {
        alpha.new = alpha.known
        mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
    }
    alpha.new = alpha.new*v
    
    
    
    
    par.old = c(  log(par$cpl[1]) , par$cpl[2])
    #a = optim( par.old, loglikgig4, ABC=ABC)
    par.ol = c(exp(par.old[1]), par.old[2])
    a = updateolU(ol=par.ol, ABC=ABC, n=2)
    #if ( loglikgig4(par.old, ABC) < loglikgig4(par.old, ABC) ) cpl.new = par$cpl
    #  else cpl.new =  c( rep(exp(a$par[1]), 2), a$par[2])
    #if(a[1]<=0){a[1]=eps}
    cpl.new =  c( a[1], a[2])
    
    
    
    
    
    
    #####second estep
    d=ncol(x)
    q=ncol(par$Lambda)
    var=par$Lambda
    fi=par$err
    fi1=ginv(fi)
    fi1=diag(d)*as.vector(fi1)
    dia=diag(q)
    if( length(q)==0){dia=q}
    inv=fi1-fi1%*%var%*%ginv(dia+t(var)%*%fi1%*%var)%*%t(var)%*%fi1
    beta=t(var)%*%inv
    term1=sweep(x,2, FUN ="-", STATS=mu.new)-abc[,1]%*%t(alpha.new)
    uhat=beta%*%t(term1)
    euu=diag(q)-beta%*%var+beta%*%t(term1)%*%term1%*%t(beta)
    
    va.new=(t(term1*(abc[,2]*weights))%*%t(uhat)%*%ginv(euu))/sum(abc[,2]*weights)
    
    A = cov.wt(x, wt=abs(abc[,2]*weights), center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data
    r = apply(x, 2, weighted.sum, wt=weights)/sumw
    R = A - (outer( r - mu.new, alpha.new) + outer(alpha.new, r - mu.new)) + outer(alpha.new,alpha.new)*ABC[1]
    
	
    fi.new=diag(R-(t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)-sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    #fi.new=diag((t(abc[,2]*term1*weights)%*%(term1)-t(2*abc[,2]*term1*weights)%*%t(uhat)%*%t(var)+sum(abc[,2]*weights)*(var%*%euu%*%t(var)))/sumw)
    #fi.new=1/fi.new
	
	
    #va.new=R%*%t(inv)%*%var
    #inv=fi1-fi1%*%va.new%*%ginv(dia+t(va.new)%*%fi1%*%va.new)%*%t(va.new)%*%fi1
    #fi.new=(R%*%t(inv))*diag(d)%*%fi
    
    var=va.new%*%t(va.new)+diag(d)*fi.new
    
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=var, cpl=cpl.new, Lambda=va.new, err=fi.new )
    
    return(new.par)
}










MAPFA <- function(data, gpar, label=NULL) {
    w = weightsFA(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}




gig2FA <- function(x=NULL, par=NULL, invS=NULL) {
    #  returns a matrix with dim length(a) x 3
    d = length(par$mu)
    
    #if (is.null(invS)) invS = ginv(par$sigma)
    alpha = par$alpha
    omega = exp(  log(par$cpl[1])  )
    
    a1 = omega + as.numeric( alpha %*% invS %*% alpha )
    b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
    v1 = par$cpl[2]-d/2
    
    val = gigFA(b=b1,a=a1, v=v1)
    return(val)
}

gigFA <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3 stima le yi
    sab =  sqrt(abs(a*b))
    kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
    kv  = besselK( sab, nu=v, expon.scaled =TRUE)
    kv12 = kv1/kv
    #sun=is.nan(kv12)
    #kv12[sun]=1
    
    sb.a = sqrt(abs(b/a))
    w    = kv12*sb.a
    invw = kv12*1/sb.a - 2*v/b
    
    #sqr1w=sqrt(sb.a)*besselK( sab, nu=v-0.5, expon.scaled =TRUE)/kv
    logw = log(sb.a) + grad( logbesselKv, x=rep(abs(v),length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    val =cbind(w,invw,logw)#,sqr1w,w2)
    return(val)
}






loglikgig4 <- function(ol=NULL, ABC=ABC) {
    val = numeric(3)
    omega  =  exp(ol[1])
    val[1] = -log(besselK( x=omega, nu=ol[2], expon.scaled=TRUE)) + omega
    val[2] = (ol[2]-1)*ABC[3]
    val[3] = -1/2* omega*(ABC[2] + ABC[1])
    val = -1*sum(val)
    return(val)
}



###Function2

##for funtion 2





ddghypFA <- function(x, par, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
    # generalized hyperbolic distribution
    
    mu    = par$mu
    sigma = par$sigma
    alpha = par$alpha
    
    #check for input error
    if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
    d = length(mu)
    omega = exp(log(par$cpl[1]))
    lambda = par$cpl[2];
    
    
    #distances between points and mu
    if (is.null(invS)) invS = ginv(sigma)
    pa = omega + alpha %*% invS %*% alpha
    mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
    
    kx = sqrt(abs(mx*pa))####################################################################################################################
    xmu = sweep(x,2,mu)
    
    lvx = matrix(0, nrow=nrow(x), 4)
    lvx[,1] = (lambda - d/2)*log(kx)
    lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
    lvx[,3] = xmu %*% invS %*% alpha
    
    
    lv = numeric(6)
    
    lv[1] = -1/2*log(abs(det( sigma ))) -d/2*(log(2)+log(pi))##################################################################################
    lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
    lv[3] = -lambda/2*( log(1) )
    lv[4] = lambda*log(omega) *0
    lv[5] = (d/2 - lambda)*log(abs( pa) )############################################################################################
    
    val = apply(lvx,1,sum) + sum(lv)
    if (!log) val = exp( val )
    
    return(val)
}



#### combinewk see Function for igpar




##end for function 2


###Function3


llikFA <- function(data,gpar, delta=0) {
    ##log likelyhood estimation
    logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) logz[,k] = ddghypFA(x=data, par=gpar[[k]], log=TRUE)
    val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
    return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                              GHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



###Function1

igpar <- function(data=NULL, g=NULL, method="kmeans") {
    ##initialization
    gpar = igpar3(data=data, G=g, n=10, method=method)
    return(gpar)
}


####Function for igpar
igpar3 <- function(data=NULL, G=NULL, n=10,label=NULL, method="kmeans") {
    #set.seed(142857)
    if(method=="modelBased"){
        l=gpcm(data,  G=G, mnames=c("VVV"))$map
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
        gpar  = rgparGH(data=data, g=G, w=z)
        
        for (j in 1:n) try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= label,  w=z)}, TRUE)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=G))
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
        gpar  = rgparGH(data=data, g=G, w=z)
        
        for (j in 1:n) try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= label,  w=z)}, TRUE)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        l=round(runif(nrow(data))*(G-1)+1)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
        gpar  = rgparGH(data=data, g=G, w=z)
        
        #l=BAR(data,l)
    }
    
    else{ l=kmeans(data,G)$cluster
        
        lw = l#kmeans(data, centers=G, iter.max=10)$cluster
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
        gpar  = rgparGH(data=data, g=G, w=z)
        
        for (j in 1:n) try({ gpar = EMgrstepGH(data=data, gpar=gpar, v=1, label= label,  w=z)}, TRUE)}
    return(gpar)
}


rgparGH <- function(data, g=2, w=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
    val = list()
    for (k in 1:g) val[[k]] = ipar(data=data, wt=w[,k])
    val$pi = rep(1/g,g)
    return(val)
}


ipar <- function(data, wt=NULL) {
    if (is.null(wt)) wt = rep(1,nrow(data))
    p = ncol(data)
    val = list()
    val$mu    = rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )
    val$alpha = rnorm(p, 0, sqrt(1/nrow(data)) )
    val$sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(val$sigma)=abs(diag(val$sigma))
    if(p==1){val$sigma=var(data)}
    
    for(i in 1:p){if(val$sigma[i,i]<0.1){val$sigma[i,i]=0.1}}
    if (any(eigen(val$sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    val$cpl   = c(1,-1/2)
    return(val)
}
#####end function for igpar

EMgrstepGH <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
    ##parameter estimtion
    if (is.null(w)) w = try(weightsGH(data=data, gpar=gpar,v=v))
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    
    G= length(gpar$pi);##number clusters
    d= length(gpar[[1]]$mu);##number variables
    Sk = array(0, c(d,d,G) )
    for (k in 1:G ) {
        #	print(gpar[[k]]$cpl)
        gpar[[k]] = updatemaScpl(x=data, par=gpar[[k]], weights=w[,k], invS=NULL, alpha.known=NULL, v=v)
        Sk[,,k]   = gpar[[k]]$sigma
        
        #		print(gpar[[k]]$cpl)
        #		print(c("******"))
    }
    gpar$pi = apply(w,2,mean)
    return(gpar)
}







updatemaScpl <- function(x, par, weights=NULL, invS=NULL, alpha.known=NULL, v=NULL) {
    ####computation of mu and alpha sigma cpl
    ##########piu importante intervenire qui!!
    if (is.null(weights)) weights=rep(1,nrow(x))
    if (is.null(invS)) invS=ginv(par$sigma)
    
    # expectations of w, 1/w, log(w) given x
    abc = gig2GH(x=x, par=par, invS=invS)
    d = length(par$mu)
    
    sumw = sum(weights)
    ABC = apply(abc,2,weighted.sum, wt=weights)/sumw
    if ( is.null(alpha.known) ) {
        A = ABC[1]
        B = ABC[2]
        u = (B - abc[,2])*weights
        t = (A*abc[,2]-1)*weights
        T = sum(t)
        
        mu.new    = apply(x, 2, weighted.sum, w=t)/T
        alpha.new = apply(x, 2, weighted.sum, w=u)/T
    } else {
        alpha.new = alpha.known
        mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
    }
    alpha.new = alpha.new*v
    
    A = cov.wt(x, wt=abc[,2]*weights, center=mu.new, method="ML")$cov*ABC[2] #returns a list containing weighted covariance matrix and the mean of the data
    #	R = A - outer(alpha.new,alpha.new)*ABC[1]
    #	r = sweep(x, 2, mu.new)
    r = apply(x, 2, weighted.sum, wt=weights)/sumw
    R = A - (outer( r - mu.new, alpha.new) + outer(alpha.new, r - mu.new)) + outer(alpha.new,alpha.new)*ABC[1]
    for(i in 1:ncol(R)){if(R[i,i]<0.00001){R[i,i]=0.00001}}
    par.old = c(  log(par$cpl[1]) , par$cpl[2])
    #a = optim( par.old, loglikgig4, ABC=ABC)
    par.ol = c(exp(par.old[1]), par.old[2])
    a = updateol(ol=par.ol, ABC=ABC, n=2)
    #if ( loglikgig4(par.old, ABC) < loglikgig4(par.old, ABC) ) cpl.new = par$cpl
    #  else cpl.new =  c( rep(exp(a$par[1]), 2), a$par[2])
    
    cpl.new =  c( a[1], a[2])
    new.par = list(mu=mu.new, alpha=alpha.new, sigma=R, cpl=cpl.new )
    return(new.par)
}






gig2GH <- function(x=NULL, par=NULL, invS=NULL) {
    #  returns a matrix with dim length(a) x 3
    d = length(par$mu)
    
    #if (is.null(invS)) invS = ginv(par$sigma)
    alpha = par$alpha
    omega = exp(  log(par$cpl[1])  )
    
    a1 = omega + as.numeric( alpha %*% invS %*% alpha )
    b1 = omega + as.numeric(mahalanobis(x, center=par$mu, cov=invS, inverted=TRUE))
    v1 = par$cpl[2]-d/2
    
    val = gigGH(b=b1,a=a1, v=v1)
    return(val)
}

gigGH <- function(a=NULL,b=NULL,v=NULL) {
    # returns a matrix with dim length(a) x 3 stima le yi
    sab =  sqrt(a*b)
    kv1 = besselK( sab, nu=v+1, expon.scaled =TRUE)
    kv  = besselK( sab, nu=v, expon.scaled =TRUE)
    kv12 = kv1/kv
    
    sb.a = sqrt(b/a)
    w    = kv12*sb.a
    invw = kv12*1/sb.a - 2*v/b
    logw = log(sb.a) + grad( logbesselKv, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    
    
    val =cbind(w,invw,logw)
    return(val)
}





###Function2


##for funtion 2


weightsGH <- function(data=NULL, gpar=NULL, v=1) {
    G = length(gpar$pi)
    if (G > 1) {
        zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
        for (k in 1:G )zlog[,k] =  ddghypGH(x=data, par=gpar[[k]], log=TRUE)
        w = t(apply(zlog, 1, function(z,wt,v) {
            fstar=v*(z + log(wt))-max(v*(z + log(wt)))
            x=exp(fstar);
            
            if (sum(x)  == 0) x= rep(1,length(x))
            x =  x/sum(x)
            return( x )
        }, wt=gpar$pi,v=v ))
    } else w = matrix(1,nrow=nrow(data), ncol=G)
    return(w)
}





ddghypGH <- function(x, par, log=FALSE, invS=NULL) {
    # x  is a n * p matrix
    # generalized hyperbolic distribution
    
    mu    = par$mu
    sigma = par$sigma
    alpha = par$alpha
    
    #check for input error
    if (length(mu) != length(alpha) ) stop("mu and alpha do not have the same length")
    d = length(mu)
    omega = exp(log(par$cpl[1]))
    lambda = par$cpl[2];
    
    
    #distances between points and mu
    if (is.null(invS)) invS = ginv(sigma)
    pa = omega + alpha %*% invS %*% alpha
    mx = omega + mahalanobis(x, center=mu, cov=invS, inverted=TRUE)
    
    kx = sqrt(mx*pa)
    xmu = sweep(x,2,mu)
    
    lvx = matrix(0, nrow=nrow(x), 4)
    lvx[,1] = (lambda - d/2)*log(kx)
    lvx[,2] = log(besselK( kx, nu=lambda-d/2, expon.scaled =TRUE)) - kx
    lvx[,3] = xmu %*% invS %*% alpha
    
    
    lv = numeric(6)
	if(is.nan(log(det( sigma )))){sigma =diag(ncol(mu))}
    
    lv[1] = -1/2*log(det( sigma )) -d/2*(log(2)+log(pi))
    lv[2] =  omega - log(besselK( omega, nu=lambda, expon.scaled =TRUE))
    lv[3] = -lambda/2*( log(1) )
    lv[4] = lambda*log(omega) *0
    lv[5] = (d/2 - lambda)*log( pa )
    
    val = apply(lvx,1,sum) + sum(lv)
    if (!log) val = exp( val )
    
    return(val)
}


##end for function 2


###Function3


llikGH <- function(data,gpar, delta=0) {
    ##log likelyhood estimation
    logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:length(gpar$pi) ) logz[,k] = ddghypGH(x=data, par=gpar[[k]], log=TRUE)
    val = sum(log(apply(logz,1,function(z,wt=NULL) {
        return(sum(exp(z)*wt))
    },wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
    return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )


###function 4

MAPGH <- function(data, gpar, label=NULL) {
    w = weightsGH(data=data, gpar=gpar, v=1)
    if (!is.null(label)) w = combinewk(weights=w, label= label)
    z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
    z = as.numeric(z)
    return( z)
}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                            MSGHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


rmgparMS <- function(g=NULL, p=NULL,data=NULL, method="kmeans",n=10) {
	val = list()
	#set.seed(142857)
    if(method=="modelBased"){
        l=gpcm(data,  G=g, mnames=c("VVV"))$map
        # l=BAR(data,l)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        #l=BAR(data,l)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        l=round(runif(nrow(data))*(g-1)+1)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=l)
    }
    else{ l=kmeans(data,g)
        l=l$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMS(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMS(data=data, gpar=gpar, v=1, label= l,  w=z)}
    return(gpar)
    return(l)
}



##########Inizialization
rgparMS <- function(g=NULL,data, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
	val = list()
    # l=rmgparMS(g, p,data,method=method)
	for (k in 1:g) val[[k]] = rparMS(data=data,wt=w[,k],k)
	val$pi = rep(1/g,g)
	return(val)
}

rparMS <- function(data,wt,k) {
	par = list()
    p=ncol(data)
	par$mu =rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )#l[k,] #rnorm(p,0,.01);
	#par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
	par$alpha =  rnorm(p, 0, sqrt(1/nrow(data)) )#rnorm(p,0,.01);
    
    sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
	par$cpl = cbind( rep(1,p), rep(-1/2,p))
    #	par$gam   = diag( rep(1,p) )
	#par$gam   = eigen( cov(matrix(rnorm((p+1)*p), p+1, p)))$vectors
    
	return(par)
}


#########Main

EMgrstepMS <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
	if (is.null(w)) w = weightsMS(data=data, gpar=gpar,v=v)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
    
	G= length(gpar$pi);
	for (k in 1:G ) {
		w = weightsMS(data=data, gpar=gpar,v=v)
		if (!is.null(label)) w = combinewk(weights=w, label= label)
		gpar[[k]] = updatemaScplpMS(y=data, par=gpar[[k]], weights=w[,k], alpha.known=NULL, v=v)
	}
	gpar$pi = apply(w,2,mean)
	return(gpar)
}



weightsMS <- function(data=NULL, gpar=NULL, v=1) {
	G = length(gpar$pi)
	if (G > 1) {
		zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
		for (k in 1:G ) zlog[,k] =  dmsghypMS(y=data, par=gpar[[k]], log=TRUE)
		w = t(apply(zlog, 1, function(z,wt,v) {
			x=exp(v*(z + log(wt)) );
            
			if (sum(x)  == 0) x= rep(1,length(x))
			x =  x/sum(x)
			return( x )
        }, wt=gpar$pi,v=v ))
	} else w = matrix(1,nrow=nrow(data), ncol=G)
	return(w)
}




dmsghypMS <- function(y, par, log=FALSE) {
    # x is a n x p matrix
    x     = y %*% (par$gam)
    mu    = par$mu;
    sigma = par$phi;
    alpha = par$alpha;
    d = length(mu); chi = par$cpl[,1]; psi = par$cpl[,1];
    lambda = par$cpl[,2];
    
    xmu = sweep(x,2,mu,"-")
    
    pa = psi + alpha^2/sigma   # p numeric
    mx = sweep(sweep(xmu^2,2,1/sigma, FUN="*"), 2,chi, "+") # n x p matrix
    kx = sqrt(sweep(mx, 2, pa, "*")) # nxp matrix
    
    lx1 = sweep( sweep(log(mx),2,log(pa),"-"), 2, (lambda - 1/2)/2, "*")
    lx2 = t(apply(kx, 1, function(z,lam=NULL) { log(besselK( z, nu=lambda-1/2, expon.scaled =TRUE)) - z }, lam=lambda ))
    lx3 = sweep(xmu, 2, alpha/sigma, FUN="*")
    
    lv = numeric(3)
    lv1 = -1/2*log( sigma ) -1/2*(log(2)+log(pi)) # d=1
    lv2 = - (log(besselK( sqrt(chi*psi), nu=lambda, expon.scaled =TRUE)) - sqrt(chi*psi) )
    lv3 = lambda/2*(log(psi)-log(chi) )
    if(ncol(y)==1){lx2=t(lx2)}
    lx=lx1 +lx2 + lx3
    #val = apply(lx1 + lx2 + lx3, 1, sum) + sum(lv1 + lv2 + lv3)
    val = apply(lx, 1, sum) + sum(lv1 + lv2 + lv3)
    if (!log) val = exp( val )
    
    return(val)
}




gig2pMS <- function(sr2=NULL, par=NULL) {
	# returns the same as gig
    
	omega = exp(log(par$cpl[,1]) )
	a1 = omega + par$alpha^2/par$phi   # vector
	B1 = sweep(sr2, 2, omega, "+") # matrix
	v1 = par$cpl[,2]-1/2 # vector
    
	val = gigp(a=a1, B=B1, v=v1)
	return(val)
}



updatemaScplpMS <- function(y=NULL, par=NULL, weights=NULL, alpha.known=NULL, v=1) {
	if (is.null(weights)) weights=rep(1,nrow(x))
	
	x = y %*% (par$gam)
	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
	abc = gig2pMS(sr2=sr2, par=par)
    
    #	new.gam= par$gam
	if (runif(1) > 1/2) new.gam = updategam2MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
	else new.gam = updategam1MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    #new.gam = update.gam1(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
	x = y %*% (new.gam)
    #	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    #print(par$cpl)
    #	abc = gig2p(sr2=sr2, par=par)
    
    
	sumw = sum(weights)
	A = apply(abc$W,   2,weighted.sum, wt=weights)/sumw
	B = apply(abc$invW,2,weighted.sum, wt=weights)/sumw
	C = apply(abc$logW,2,weighted.sum, wt=weights)/sumw
	if ( is.null(alpha.known) ) {
		u = sweep(sweep(-abc$invW, 2, B, "+"), 1, weights, "*")
		t = sweep(sweep(abc$invW, 2, A, "*")-1, 1, weights, "*")
		T = apply(t,2,sum)
        
		mu.new    = apply(x*t,2,sum)/T
		alpha.new = apply(x*u,2,sum)/T
        
	} else {
		alpha.new = alpha.known
		mu.new    = apply(x*abc$invW, 2, weighted.sum, wt=weights)/sumw - alpha.new/B
	}
	alpha.new = alpha.new*v
    
	# update for sigma
	Ax = apply(sweep(x, 2, mu.new, "-")^2*abc$invW, 2, weighted.sum, wt=weights)/sumw
	ax = apply(x, 2, weighted.sum, wt=weights)/sumw
	sigma.new = Ax - 2*(ax - mu.new)*alpha.new + alpha.new^2*A
    
	omega=  exp(log(par$cpl[,1]) )
	test  = cbind(omega, lambda=par$cpl[,2], A, B, C)
	cpl.new = t(apply(test, 1, function(z) {
		temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
		return( c( ( temp[1]), temp[2]) )
    }))
    
    #	new.gam = update.gam1(gam0=par$gam, x=x, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=par$gam, y=y, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=new.gam, x=y %*% (new.gam), sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    
	new.par = list(mu=mu.new, alpha=alpha.new,phi=sigma.new, cpl=cpl.new, gam= new.gam )
	return(new.par)
}


updategam2MS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
	x = y %*% (gam0)
	sumw = sum(wt)
	Bs = sweep(invW, 2, 1/sigma, "*" )
    #	u  = sweep(y * invW, 1, wt, "*")
	wty = sweep(y, 1, wt, "*")
	F0t = t( x * Bs ) %*% wty
    
	e2  = apply(y^2, 1, sum)*wt
    
    if(ncol(y)==1){F1t = (apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)}
    else{
        F1t = diag(apply(Bs, 2, weighted.sum, wt = e2)) %*% t(gam0)}
    
	v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
	A0 = t(v)  %*% wty
	
	F2 = F0t - F1t - A0
	z  = svd( F2/sumw)
	gam = ( (z$v) %*% t(z$u) )
	gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
	
	if (objgamMS(gam0=gam,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) <  objgamMS(gam0=gam0,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) ) {
        #		print('not good enough')
		gam = gam0
    }
	
	return(gam)
}






updategam1MS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
    
	x = y %*% (gam0)
	sumw = sum(wt)
	Bs = sweep(invW, 2, 1/sigma, "*" )
	u  = sweep(y * invW, 1, wt, "*")
	wty = sweep(y, 1, wt, "*")
	F0t = t( x * Bs ) %*% wty
    
	e2  = apply(Bs, 1, max)*wt
	F1t = (cov.wt(y, center=rep(0,ncol(y)), wt=e2, method="ML" )$cov*sum(e2)) %*% (gam0)
    
	v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
	A0 = t(v)  %*% wty
	
	F2 = F0t - t(F1t) - A0
	z  = svd( F2/sumw)
	gam = ( (z$v) %*% t(z$u) )
	gam = sweep(gam, 2, sign(diag(gam)), FUN="*" )
	
    if (objgamMS(gam0=gam,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) <  objgamMS(gam0=gam0,y=y, sigma=sigma, alpha=alpha, mu=mu, wt=wt, invW=invW) ) {
        #		print('not good enough')
		gam = gam0
    }
	
	return(gam)
}




objgamMS <- function(gam0=NULL, y=NULL, sigma=NULL, alpha=NULL, mu=NULL, wt=NULL, invW=NULL) {
	
	x = y %*% (gam0)
	sumw = sum(wt)
	Bs  = sweep(invW, 2, 1/sigma, "*" )
	wtx = sweep(x, 1, wt, "*")
	F0t = t( x * Bs ) %*% wtx
    
	v  = sweep(sweep(Bs, 2, mu, "*"), 2, alpha/sigma, "+")
	A0 = t(v)  %*% wtx
    
	val = -1*(tr(F0t ) - 2*tr( A0  ))/sumw
	
	return(val)
}


weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )


################likelihood

llikMS <- function(data, gpar) {
	logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
	for (k in 1:length(gpar$pi) ) logz[,k] = dmsghypMS(y=data, par=gpar[[k]], log=TRUE)
	val = sum(log(apply(logz,1,function(z,wt=NULL) {
		return(sum(exp(z)*wt))
	},wt=gpar$pi)))
    #	val = sum(log(apply(z,1,weighted.sum,wt=gpar$pi)))
    #	if (is.nan(val)) {
    #		print(gpar)
    #		print(logz)
    #		}
	return(val)
}


#############MAP

MAPMS <- function(data, gpar, label=NULL) {
	w = weightsMS(data=data, gpar=gpar, v=1)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
	z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
	z = as.numeric(z)
	return( z)
}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
###                                                            cMSGHD                                                          ###
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
rmgparMSr <- function(g=NULL, data=NULL, method="kmeans",n=10) {
	val = list()
	#set.seed(142857)
    if(method=="modelBased"){
        l=Mclust(data,G=g)$classification
        # l=BAR(data,l)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z)
    }
    else if(method=="hierarchical"){
        l=(cutree(hclust(dist(data),"ward.D"), k=g))
        #l=BAR(data,l)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z)
    }
    else if(method=="random"){
        #l=kmeans(data,g)
        l=round(runif(nrow(data))*(g-1)+1)
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=l)}
    else{ l=kmeans(data,g)
        l=l$cluster#$centers}
        z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=g), label=l)
        gpar  = rgparMSr(data=data, g=g, w=z,l=l)
        
        for (j in 1:n)  gpar = EMgrstepMSr(data=data, gpar=gpar, v=1, label= l,  w=z)}
    return(gpar)
    return(l)
}




##########Inizialization
####################################change
rgparMSr<- function(g=NULL,data, w=NULL,l=NULL) {
    if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
	val = list()
    # l=rmgparMS(g, p,data,method=method)
	for (k in 1:g) val[[k]] = rparMSr(data=data,wt=w[,k],k)
	val$pi = rep(1/g,g)
	return(val)
}

####################################change
rparMSr <- function(data,wt,k) {
    par = list()
    p=ncol(data)
	par$mu =rnorm(p, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )#l[k,] #rnorm(p,0,.01);
	#par$phi  =  rnorm(p, 0, sqrt(1/nrow(data)) )#rep(1,p);
	par$alpha =  rnorm(p, 0, sqrt(1/nrow(data)) )#rnorm(p,0,.01);
    
    sigma = diag( diag( cov.wt(data, wt = wt, method="ML")$cov) ) #+ diag(apply(data,2,var))*1 # + outer(val$alpha,val$alpha)
    diag(sigma)=abs(diag(sigma))
    if(p==1){sigma=var(data)}
    
    for(i in 1:p){if(sigma[i,i]<0.1){sigma[i,i]=0.1}}
    if (any(eigen(sigma)$values <= 0 ) ) val$sigma =  diag(apply(data,2,var))
    par$gam   = eigen( sigma)$vectors
    par$phi  = eigen( sigma)$values
    par$cpl = cbind( rep(1,p), rep(1,p))
	return(par)
}


#########Main

EMgrstepMSr <- function(data=NULL, gpar=NULL, v=1, label=NULL, w=NULL) {
	if (is.null(w)) w = weightsMS(data=data, gpar=gpar,v=v)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
    
	G= length(gpar$pi);
	for (k in 1:G ) {
		w = weightsMS(data=data, gpar=gpar,v=v)
		if (!is.null(label)) w = combinewk(weights=w, label= label)
        
		gpar[[k]] = updatemaScplpMSr(y=data, par=gpar[[k]], weights=w[,k], alpha.known=NULL, v=v)
	}
	gpar$pi = apply(w,2,mean)
	return(gpar)
}


####################################change
updatemaScplpMSr <- function(y=NULL, par=NULL, weights=NULL, alpha.known=NULL, v=1) {
	if (is.null(weights)) weights=rep(1,nrow(x))
	
	x = y %*% (par$gam)
	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
	abc = gig2pMS(sr2=sr2, par=par)
    
    #	new.gam= par$gam
	if (runif(1) > 1/2) new.gam = updategam2MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
	else new.gam = updategam1MS(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
    #new.gam = update.gam1(gam0=par$gam, y=y, sigma=par$phi, alpha=par$alpha, mu=par$mu, wt=weights, invW= abc$invW)
	x = y %*% (new.gam)
    #	sr2= sweep( sweep(x,2,par$mu,"-")^2, 2, 1/par$phi, FUN="*")
    #print(par$cpl)
    #	abc = gig2p(sr2=sr2, par=par)
    
    
	sumw = sum(weights)
	A = apply(abc$W,   2,weighted.sum, wt=weights)/sumw
	B = apply(abc$invW,2,weighted.sum, wt=weights)/sumw
	C = apply(abc$logW,2,weighted.sum, wt=weights)/sumw
	if ( is.null(alpha.known) ) {
		u = sweep(sweep(-abc$invW, 2, B, "+"), 1, weights, "*")
		t = sweep(sweep(abc$invW, 2, A, "*")-1, 1, weights, "*")
		T = apply(t,2,sum)
        
		mu.new    = apply(x*t,2,sum)/T
		alpha.new = apply(x*u,2,sum)/T
        
	} else {
		alpha.new = alpha.known
		mu.new    = apply(x*abc$invW, 2, weighted.sum, wt=weights)/sumw - alpha.new/B
	}
	alpha.new = alpha.new*v
    
	# update for sigma
	Ax = apply(sweep(x, 2, mu.new, "-")^2*abc$invW, 2, weighted.sum, wt=weights)/sumw
	ax = apply(x, 2, weighted.sum, wt=weights)/sumw
	sigma.new = Ax - 2*(ax - mu.new)*alpha.new + alpha.new^2*A
    
	omega=  exp(log(par$cpl[,1]) )
	test  = cbind(omega, lambda=par$cpl[,2], A, B, C)
	cpl.new = t(apply(test, 1, function(z) {
		temp = updateol(ol= z[1:2], ABC=z[3:5], n=2)
		return( c( ( temp[1]), temp[2]) )
    }))
    
    #	new.gam = update.gam1(gam0=par$gam, x=x, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=par$gam, y=y, sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    #	new.gam = update.gam2(gam0=new.gam, x=y %*% (new.gam), sigma=sigma.new, alpha=alpha.new, mu=mu.new, wt=weights, invW= abc$invW)
    for(i in 1:nrow(cpl.new)){
        if( cpl.new[i,2]<1){cpl.new[i,2]=1}}
	new.par = list(mu=mu.new, alpha=alpha.new, phi=sigma.new, cpl=cpl.new, gam= new.gam )
	return(new.par)
}
