EthanolDiff_BTCS<-function(Vol,U_thresh,eps){
  #Ethanol diffusion in tumor
  # Vol = 0.2
  # U_thresh = 0.25
  # eps=1
  dr = 0.01
  dt = dr
  nu = dt/dr
  mu = dt/dr^2
  T = 2
  X = 1
  N = round(T/dt + 1)
  M = round(X/dr + 1)
  
  M_kill = matrix(0,N,1)
  Percent_kill = 0
  r = seq(0,1,by = dr)
  
  U = matrix(0,M,N)
  A = matrix(0,M,M)
  
  width = 1
  r1 = round(length(r)*width)
  height = Vol/(sum((r[1:(r1 - 1)]^2)*exp(1 - 1./(1 - (r[1:(r1 - 1)]/r[r1])^2)))*dr*4*pi)
  U[1:r1 - 1,1] = exp(1 - 1/(1 - (r[1:(r1 - 1)]/r[r1])^2))*(height)
  #U(:,1) = height*(1-r/2);

  for(n in 2:N){
    #Point at origin
    A[1,1] = mu + 1
    A[1,2] = -mu
    A[M,M] = mu + 1 + eps*dr*(nu + mu)
    A[M,M-1] = -mu
    
    #Points in the middle
    
    for (m in  2:(M-1)){
      A[m,m+1] = -nu/(2*r[m]) - mu/2
      A[m,m] = mu + 1
      A[m,m-1] = nu/(2*r[m]) - mu/2
      
      if (U[m,n-1] >= U_thresh){
        M_kill[n-1] = r[m]^3
      }
    }
      
    U[,n] = solve(A)%*%U[,n-1]
      
    if (U[M,n-1] >= U_thresh){
      M_kill[n-1] = 1
      break
    }
  }
  return(max(M_kill))
}
