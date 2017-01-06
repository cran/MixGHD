
contourpl<-function(input){
  z=input@data
  if(input@scale==T){
  z=scale(z)}
  method=input@method
  xmin=min(z[,1])-1
  xmax=max(z[,1])+1
  ymin=min(z[,2])-1
  ymax=max(z[,2])+1
  em=input@gpar
  G=length(em$pi)
  x = seq(xmin,xmax,length.out=50)
  y = seq(ymin,ymax,length.out=50)
  
  
  xyS1 = array(0,c(length(x),length(y),G))
  if(method=="MCGHD"){
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  dmsghyp(xy, em[[k]], log=FALSE) 
          
        }
      }
    }
  }
  else if(method=="MGHD"){
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  ddghypGH(xy, em[[k]], log=FALSE) 
          
        }
      }
    }
  }
  else{
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        for(k in 1:G){
          xy <- matrix(cbind(x[i],y[j]),1,2)	
          xyS1[i,j,k] =  dmsghypMS(xy, em[[k]], log=FALSE) 
          
        }
      }
    }
  }
  zz=xyS1[,,1]*em$pi[1]
  for(k in 2:G){
    zz=xyS1[,,k]*em$pi[k]+zz}
  
  contour(x=x,y=y,z=zz, levels=c(.005,.01,.025,.05, .1,.25), main=method,ylim=c(ymin,ymax), xlim=c(xmin,xmax))
  points(z,col=input@map,pch=input@map)}