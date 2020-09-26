hourglass = function(kh1=0.01, kh2=0.007, kh3=0.012, v=0.014, ini.c=1, ini.box=1,ntimesteps=200){
  
  #transition matrix:K0
  K = read.csv("hourglass.csv")
  for(i in 1:50){
    for(j in 1:50){
      if(K[i,j]=="-2*kh1-v" ) {K[i,j]=-2*kh1-v}
      else if(K[i,j]=="kh1" ) {K[i,j]=kh1}
      else if(K[i,j]=="kh2" ) {K[i,j]=kh2}
      else if(K[i,j]=="kh3" ) {K[i,j]=kh3}
      else if(K[i,j]=="2*kh2" ) {K[i,j]=2*kh2}
      else if(K[i,j]=="2*kh3" ) {K[i,j]=2*kh3}
      else if(K[i,j]=="kh1+v" ) {K[i,j]=kh1+v}
      else if(K[i,j]=="-2*kh1-kh2-v" ) {K[i,j]=-2*kh1-kh2-v}
      else if(K[i,j]=="-2*kh1-kh3-v" ) {K[i,j]=-2*kh1-kh3-v}
      else if(K[i,j]=="-2*kh1-2*kh2-v" ) {K[i,j]=-2*kh1-2*kh2-v}
      else if(K[i,j]=="-2*kh1-2*kh3-v" ) {K[i,j]=-2*kh1-2*kh3-v}
      else if(K[i,j]=="-2*kh1-kh2-kh3-v" ) {K[i,j]=2*kh1-kh2-kh3-v}
      else if(K[i,j]=="-2*kh2-4*kh3" ) {K[i,j]=-2*kh2-4*kh3}
      else if(K[i,j]=="-4*kh2-2*kh3" ) {K[i,j]=-4*kh2-2*kh3}
      else K[i,j]=0
    }
  }
 as.matrix(K)
 K0=matrix(as.numeric(unlist(K)),nrow=nrow(K))
  
 #find eigen
  eigen = eigen(K0)
  Xe = eigen$vectors
  lamda = eigen$values 
  
}
  
  
