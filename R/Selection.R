Selection<-function(data,label=NULL, G=2, niter=50,method="GHD",starting="kmeans",max.iter=10,eps=1e-2,labelcl=NULL,q=2,scale=TRUE,criterion="ARI"){
   if( criterion =="ARI"){
       if (is.null(label)) stop('label is null, please use criterion BIC ')
  ad=0
  mod=0
  if(method=="MGHFA"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MGHFA(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
			if(is.list(res)){it=it+1}}
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
      
    }}
  else if(method=="MSGHD"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
			if(is.list(res)){it=it+1}}
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else if(method=="cMSGHD"){
      for(i in 1:niter){
		  it=niter-1
		  while(it<niter){
          res=try(cMSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
		  if(is.list(res)){it=it+1}}
          adn=adjustedRandIndex(res$map,label)
          if(adn>ad){
              mod=res
              ad=adn}
      }}
  else if(method=="MCGHD"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
			 if(is.list(res)){it=it+1}}
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else {
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
		 if(is.list(res)){it=it+1}}
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}}
   ### ICL
   else if ( criterion =="ICL"){
       mod=0
       ICL=-Inf
       if(method=="MGHFA"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   res=try(MGHFA(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               icln=res$ICL
               if(icln>ICL){
                   mod=res
                   ICL=icln}
               
           }}
       else if(method=="MSGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   res=try(MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       else if(method=="cMSGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   res=try(cMSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       else if(method=="MCGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   
                   res=try(MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}}}
       else {
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   
                   res=try(MGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       if(is.null(label)){ad=NA}
       else{
           ad=adjustedRandIndex(res$map,label)}
       
       
   }

   #######BIC
    else{
         mod=0
        BIC=-Inf
        if(method=="MGHFA"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(MGHFA(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
               if(is.list(res)){it=it+1}}
					bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
                
            }}
        else if(method=="MSGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                if(is.list(res)){it=it+1}}
				bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="cMSGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(cMSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
					 if(is.list(res)){it=it+1}}
                bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="MCGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){

                res=try(MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                if(is.list(res)){it=it+1}}
				bicn=res$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}}}
        else {
            for(i in 1:niter){
				it=niter-1
				while(it<niter){

                res=try(MGHD(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                 if(is.list(res)){it=it+1}}
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