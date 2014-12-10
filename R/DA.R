DA<-function(train,trainL,test,testL,method="GHD",starting="kmeans",max.iter=10,q=2,scale=TRUE){
    if (min(trainL)==0) stop('training label equal to 0')

    G=max(trainL)
  if(method=="MGHFA"){

      model=MGHFA(train,G=G,method=starting,label=trainL,q=q,scale=scale)
      testmodel=MAPFA(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,trainL)
      aritest=adjustedRandIndex(testmodel,testL)
      
    }
  else if(method=="MSGHD"){

    model=MSGHD(train,G=G,max.iter=max.iter,method=starting,label=trainL,scale=scale)
      testmodel=MAPMS(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,trainL)
      aritest=adjustedRandIndex(testmodel,testL)  }
  else if(method=="MCGHD"){
   
      model=MCGHD(train,G=G,max.iter=max.iter,method=starting,label=trainL,scale=scale)
      testmodel=MAP(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,trainL)
      aritest=adjustedRandIndex(testmodel,testL)
      }
  else if(method=="cMSGHD"){
      
      model=cMSGHD(train,G=G,max.iter=max.iter,method=starting,label=trainL,scale=scale)
      testmodel=MAPMS(scale(test),gpar=model$gpar)
      aritrain=adjustedRandIndex(model$map,trainL)
      aritest=adjustedRandIndex(testmodel,testL)
  }
  else {

      res=MGHD(train,G=G,n=max.iter,method=starting,label=trainL,scale=scale)
      testmodel=MAPGH(scale(test),gpar=res$gpar)
      aritrain=adjustedRandIndex(res$map,trainL)
      aritest=adjustedRandIndex(testmodel,testL)
  }
  return(list(model=model,testMembership=testmodel,ARItest=aritest,ARItrain=aritrain))
}