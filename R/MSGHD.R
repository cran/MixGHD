
MSGHD <- function(data=NULL, gpar0=NULL, G=2, max.iter=100, label =NULL ,eps=1e-2, method="km",scale=TRUE ) {
    data=as.matrix(data)
    if( scale==TRUE)
    {data=scale(data)}
    pcol=ncol(data)
    #  if (nrow(data)<((G-1)+G*(4*pcol+pcol*(pcol-1)/2)))stop('G is too big, number of parameters > n')
    if (is.null(data)) stop('data is null')
    if (nrow(data) == 1) stop('nrow(data) is equal to 1')
    # if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
    if (any(is.na(data))) stop('No NAs allowed.')
    if (is.null(G)) stop('G is NULL')
    # if ( G < 1) stop('G is not a positive integer')
    if (  max.iter< 1) stop('max.iter is not a positive integer')
    
    
    
    bico=-Inf
    t=length(G)
    BIC=matrix(NA,t,1)
    cont=0
	for(b in 1:t){
        mo=try(mainMSGHD(data=data, gpar0=gpar0, G=G[b], n=max.iter, eps=eps,  label=label,method= method),silent = TRUE)
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
    val=list(BIC=BIC,model=model)
        cat("The best model (BIC) for the range of  components used is  G = ", sg,".\nThe BIC for this model is ", bico,".",sep="")
        return(val)
}
