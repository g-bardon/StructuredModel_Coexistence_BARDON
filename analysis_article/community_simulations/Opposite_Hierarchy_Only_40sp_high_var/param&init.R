library("mvtnorm")
init <- function(S){
  nj <- runif(S, min=10, max=50)
  na <- runif(S, min=10, max=50)
  return(cbind(nj,na))
}

parametrage <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.5)   ## log-normal distribution for fertility, low variance to avoid too large fertility
  mat_rate = rbeta(S, 5,5) ## beta distribution centered in 1/2
  juv_surv = rbeta(S, 5,5) ## beta distribution centered in 1/2
  adult_surv = rbeta(S, 5,5) ## beta distribution centered in 1/2
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Normal distribution for competition coefficients 
  beta_intra = rnorm(S, 0.05, 0.01)
  sig=0.02
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.04,0.04), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)  ## Matrix of alpha competition coefficients : i column corresponds to the alpha competition undergoes by species i
  beta=matrix(0,S,S)  ## Matrix of beta competition coefficients : i column corresponds to the beta competition undergoes by species i
  
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## intraspecific competition coefficients on the diagonal
        beta[i,j] <- beta_intra[i]
      }else{
        x <- alpha_beta_inter[c,] 
        alpha[i,j] <- x[1] ## alpha and beta interspecific competition coefficient with negative correlation
        beta[i,j] <- x[2]
        c=c+1
      }
    }
  }
  
  return(list(fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Return the parameter set with the negative correlation structure
}

S=40
param = parametrage(S)
ini <- init(S)

df_param_vital_rates <- data.frame("pi"=param[[1]], "gamma"=param[[2]], "phi"=param[[3]], "sa"=param[[4]])
write.csv(df_param_vital_rates, file = "param_vital_rates.csv", row.names = FALSE)
write.table(param[[5]], file="matrix_alpha.csv", row.names = FALSE, col.names = FALSE)
write.table(param[[6]], file="matrix_beta.csv", row.names = FALSE, col.names = FALSE)
write.csv(ini, file = "initialization.csv", row.names = FALSE)
