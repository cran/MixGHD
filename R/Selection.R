
Selection<-function(data,label=NULL, G=2, niter=50,method="GHD",starting="kmeans",max.iter=10,eps=1e-2,labelcl=NULL,q=2,scale=TRUE,criterion="ARI"){
   if( criterion =="ARI"){
       if (is.null(label)) stop('label is null, please use criterion BIC ')
  ad=0
  mod=0
  if(method=="MGHFA"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MGHFASel(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
			if(is.list(res)){it=it+1}}
      adn=ARI(res$model$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
      
    }}
  else if(method=="MSGHD"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
			if(is.list(res)){it=it+1}}
      adn=ARI(res$model$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else if(method=="cMSGHD"){
      for(i in 1:niter){
		  it=niter-1
		  while(it<niter){
          res=try(cMSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
		  if(is.list(res)){it=it+1}}
          adn=ARI(res$model$map,label)
          if(adn>ad){
              mod=res
              ad=adn}
      }}
  else if(method=="MCGHD"){
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MCGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
			 if(is.list(res)){it=it+1}}
      adn=ARI(res$model$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else {
    for(i in 1:niter){
		it=niter-1
		while(it<niter){
      res=try(MGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
		 if(is.list(res)){it=it+1}}
      adn=ARI(res$model$map,label)
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
                   res=try(MGHFASel(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               icln=res$model$ICL
               if(icln>ICL){
                   mod=res
                   ICL=icln}
               
           }}
       else if(method=="MSGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   res=try(MSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$model$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       else if(method=="cMSGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   res=try(cMSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$model$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       else if(method=="MCGHD"){
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   
                   res=try(MCGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$model$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}}}
       else {
           for(i in 1:niter){
               it=niter-1
               while(it<niter){
                   
                   res=try(MGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                   if(is.list(res)){it=it+1}}
               ICLn=res$model$ICL
               if(ICLn>ICL){
                   mod=res
                   ICL=ICLn}
           }}
       if(is.null(label)){ad=NA}
       else{
           ad=ARI(res$model$map,label)}
       
       
   }

   #######BIC
    else{
         mod=0
        BIC=-Inf
        if(method=="MGHFA"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(MGHFASel(data,G=G,method=starting,label=labelcl,max.iter=max.iter,q=q,scale=scale,eps=eps))
               if(is.list(res)){it=it+1}}
					bicn=res$model$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
                
            }}
        else if(method=="MSGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(MSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                if(is.list(res)){it=it+1}}
				bicn=res$model$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="cMSGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){
                res=try(cMSGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
					 if(is.list(res)){it=it+1}}
                bicn=res$model$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        else if(method=="MCGHD"){
            for(i in 1:niter){
				it=niter-1
				while(it<niter){

                res=try(MCGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                if(is.list(res)){it=it+1}}
				bicn=res$model$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}}}
        else {
            for(i in 1:niter){
				it=niter-1
				while(it<niter){

                res=try(MGHDSel(data,G=G,max.iter=max.iter,method=starting,label=labelcl,scale=scale,eps=eps))
                 if(is.list(res)){it=it+1}}
					bicn=res$model$BIC
                if(bicn>BIC){
                    mod=res
                    BIC=bicn}
            }}
        if(is.null(label)){ad=NA}
        else{
            ad=ARI(res$model$map,label)}
    
    
    }
    if(is.numeric(mod)){mod=res
    ad=adn}
return(l=list(model=mod,ARI=ad))
}
