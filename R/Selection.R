Selection<-function(data,label, G=2, niter=50,method="GHD",starting="kmeans",max.iter=10,labelclassification=NULL,q=2){
  ad=0
  if(method=="MGHFA"){
    for(i in 1:niter){
      res=MGHFA(data,G=G,method=starting,label=labelclassification,q=q)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
      
    }}
  else if(method=="MSGHD"){
    for(i in 1:niter){
      res=MSGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else if(method=="MCGHD"){
    for(i in 1:niter){
      res=MCGHD(data,G=G,max.iter=max.iter,method=starting,label=labelclassification)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  else {
    for(i in 1:niter){
      res=MGHD(data,G=G,n=max.iter,method=starting,label=labelclassification)
      adn=adjustedRandIndex(res$map,label)
      if(adn>ad){
        mod=res
        ad=adn}
    }}
  return(l=list(model=mod,ARI=ad))
}