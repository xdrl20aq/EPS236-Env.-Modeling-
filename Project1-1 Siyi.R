string.of.boxes = function(nbox=5, k.h=0.01, va=0.014, ini.c=1, ini.box=1,ntimesteps=200){
  
  #transition matrix:K
  K = matrix(0, ncol=nbox, nrow=nbox)
  K[1,] = c(-2*k.h-va, k.h, rep(0,nbox-3), k.h+va)
  for(n in 2:(nbox-1)) {
    K[n,] = c(rep(0,n-2), k.h+va, -2*k.h-va, k.h ,rep(0,nbox-n-1))
  }
  K[nbox,] = c(k.h, rep(0,nbox-3), k.h+va, -2*k.h-va)
  
  #find eigen
  eigen = eigen(K)
  Xe = eigen$vectors
  lamda = eigen$values 
  

  #solve equation
  c0 = c(rep(0,nbox));  c0[ini.box] = ini.c
  c0 = as.matrix(c0)
  c00 = solve(Xe, c0)
  ct = c0
  for(n in 1:ntimesteps){
    ctcol = plot(1:101, ct[1,],col = "red") diag(exp(lamda*n)) %*% c00
    cc=Xe%*%ctcol 
    ct=cbind(ct,cc)
  }

  return(ct)
}

####matplot
ct = string.of.boxes(ntimesteps=200)
plot(1:(ntimesteps+1), ct[1,],col="1",ylim=0:1)
for (i in 2:5){
  points(1:(ntimesteps+1), ct[i,],col = i)
}




