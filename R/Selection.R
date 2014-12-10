Selection<-function(data,label=NULL, G=2, niter=50,method="GHD",starting="kmeans",max.iter=10,labelclassification=NULL,q=2,scale=TRUE,criterion="ARI"){
   if( criterion =="ARI"){
       if (is.null(label)) stop('label is null, please use criterion BIC ')
  ad=0
  mod=0
  if(method=="MGHFA"){
    for(i in 1:niter){
      res=MGHFA(data,G=G,method=starting,label=labelclassification,q=q,scale=scale)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
      
    }}
  else if(method=="MSGHD"){
    for(i in 1:niter){
      res=MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else if(method=="cMSGHD"){
      for(i in 1:niter){
          res=cMSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
          adn=adjustedRandIndex(res$map,label)
          if(adn>ad){
              mod=res
              ad=adn}
      }}
  else if(method=="MCGHD"){
    for(i in 1:niter){
      res=MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else {
    for(i in 1:niter){
      res=MGHD(data,G=G,n=max.iter,method=starting,label=labelclassification,scale=scale)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}}
    else{
         mod=0
        BIC=-Inf
        if(method=="MGHFA"){
            for(i in 1:niter){
                res=MGHFA(data,G=G,method=starting,label=labelclassification,q=q,scale=scale)
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
                
            }}
        else if(method=="MSGHD"){
            for(i in 1:niter){
                res=MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="cMSGHD"){
            for(i in 1:niter){
                res=cMSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="MCGHD"){
            for(i in 1:niter){
                res=MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification,scale=scale)
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}}}
        else {
            for(i in 1:niter){
                res=MGHD(data,G=G,n=max.iter,method=starting,label=labelclassification,scale=scale)
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        if(is.null(label)){ad=NA}
        else{
            ad=adjustedRandIndex(res$map,label)}
    
    
    }
    if(is.numeric(mod)){mod=res
    ad=adn}
return(l=list(model=mod,ARI=ad))
}