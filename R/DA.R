DA<-function(training,traininglabel,test,testlabel,method="GHD",starting="kmeans",max.iter=10,q=2){
    if (min(traininglabel)==0) stop('training label equal to 0')

    G=max(traininglabel)
  if(method=="MGHFA"){

      model=MGHFA(training,G=G,method=starting,label=traininglabel,q=q)
      testmodel=MAPFA(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,traininglabel)
      aritest=adjustedRandIndex(testmodel,testlabel)
      
    }
  else if(method=="MSGHD"){

    model=MSGHD(training,G=G,max.iter=max.iter,method=starting,label=traininglabel)
      testmodel=MAPMS(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,traininglabel)
      aritest=adjustedRandIndex(testmodel,testlabel)  }
  else if(method=="MCGHD"){
   
      model=MCGHD(training,G=G,max.iter=max.iter,method=starting,label=traininglabel)
      testmodel=MAP(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,traininglabel)
      aritest=adjustedRandIndex(testmodel,testlabel)
      }
  else if(method=="cMSGHD"){
      
      model=cMSGHD(training,G=G,max.iter=max.iter,method=starting,label=traininglabel)
      testmodel=MAPMS(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,traininglabel)
      aritest=adjustedRandIndex(testmodel,testlabel)
  }
  else {

      res=MGHD(training,G=G,n=max.iter,method=starting,label=traininglabel)
      testmodel=MAPGH(scale(test),gpar=res$gpar)
      aritrain=adjustedRandIndex(res$map,traininglabel)
      aritest=adjustedRandIndex(testmodel,testlabel)
  }
  return(list(model=model,testMembership=testmodel,ARItest=aritest,ARItrain=aritrain))
}