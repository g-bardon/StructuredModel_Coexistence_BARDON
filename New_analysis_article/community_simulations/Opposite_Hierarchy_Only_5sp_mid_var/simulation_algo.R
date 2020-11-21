algo_simul <- function(gen, S, ini, fert, mat_rate, juv_surv, adult_surv, alpha, beta){  
  nj = ini[,1]
  na = ini[,2]
  nj_plus = rep(0, S)
  na_plus = rep(0, S)
  Rs=rep(0, S)
  Rf=rep(0, S)
  A=list()
  for (c in 1:gen){
    for (i in 1:S){
      Rs[i] = alpha[,i] %*% na 
      Rf[i] = beta[,i] %*% na 
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)
      n =A[[i]]%*% c(nj[i], na[i])
      nj_plus[i]=n[1]
      na_plus[i]=n[2]
    }
    nj=nj_plus
    na=na_plus
  }
  return(sum(na+nj>1))
}    

## Simulation with complete trajectory
algo_simul_traj <- function(gen, S, ini, fert, mat_rate, juv_surv, adult_surv, alpha, beta){  
  nj = ini[,1]
  na = ini[,2]
  nj_plus = rep(0, S)
  na_plus = rep(0, S)
  traj<- matrix(0, nrow=40, ncol=gen)
  traj[,1] <- nj+na
  Rs=rep(0, S)
  Rf=rep(0, S)
  A=list()
  for (c in 2:gen){
    for (i in 1:S){
      Rs[i] = alpha[,i] %*% na 
      Rf[i] = beta[,i] %*% na 
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)
      n =A[[i]]%*% c(nj[i], na[i])
      nj_plus[i]=n[1]
      na_plus[i]=n[2]
    }
    nj=nj_plus
    na=na_plus
    traj[,c] <- nj+na
  }
  return(traj)
} 

