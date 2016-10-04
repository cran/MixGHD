
contourpl<-function(data,output,method="MGHD",scale=TRUE){
 z=data
   if(scale){
  z=scale(z)}
  
  xmin=min(z[,1])-1
  xmax=max(z[,1])+1
  ymin=min(z[,2])-1
  ymax=max(z[,2])+1
  em=output$model
  G=length(em$gpar$pi)
  x = seq(xmin,xmax,length.out=50)
  y = seq(ymin,ymax,length.out=50)
  
  
  xyS1 = array(0,c(length(x),length(y),G))
  if(method=="MCGHD"){
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  dmsghyp(xy, em$gpar[[k]], log=FALSE) 
          
        }
      }
    }
  }
  else if(method=="MGHD"){
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  ddghypGH(xy, em$gpar[[k]], log=FALSE) 
          
        }
      }
    }
  }
  else{
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  dmsghypMS(xy, em$gpar[[k]], log=FALSE) 
          
        }
      }
    }
  }
  zz=xyS1[,,1]*em$gpar$pi[1]
  for(k in 2:G){
    zz=xyS1[,,k]*em$gpar$pi[k]+zz}
  
  contour(x=x,y=y,z=zz, levels=c(.005,.01,.025,.05, .1,.25), main=method,ylim=c(ymin,ymax), xlim=c(xmin,xmax))
  points(z,col=em$map,pch=em$map)}