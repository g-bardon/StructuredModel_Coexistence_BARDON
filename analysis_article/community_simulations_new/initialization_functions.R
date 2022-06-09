library("mvtnorm")
init <- function(S){
  nj <- runif(S, min=10, max=50)
  na <- runif(S, min=10, max=50)
  return(cbind(nj,na))
}

parametrage <- function (S, mu_fert, sig_fert, ab, mu_alpha, sig_alpha, mu_beta, sig_beta, sig_inter, mu_inter){
  fert = rlnorm(S, meanlog = mu_fert, sdlog=sig_fert)  ## log-normal distribution for fertility, low variance to avoid too large fertility
  mat_rate = rbeta(S, ab,ab) ## beta distribution centered in 1/2
  juv_surv = rbeta(S, ab,ab) ## beta distribution centered in 1/2
  adult_surv = rbeta(S, ab,ab) ## beta distribution centered in 1/2
  
  alpha_intra = rnorm(S, mu_alpha, sig_alpha)  # Normal distribution for competition coefficients 
  beta_intra = rnorm(S, mu_beta, sig_beta)
  sig=sig_inter
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(mu_inter,mu_inter), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrix of alpha competition coefficients : i column corresponds to the alpha competition undergoes by species i
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
