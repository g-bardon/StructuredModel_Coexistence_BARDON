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
      Rs[i] = alpha[,i] %*% nj
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

simul_permut <- function(gen, S, ini, fert, mat_rate, juv_surv, adult_surv, alpha, beta){
  x <- c()
  for (k in 1:200){  
    liste_beta_permu <- sample(beta[which(lower.tri(beta) | upper.tri(beta))], length(beta[which(lower.tri(beta) | upper.tri(beta))]))  ## Recup tous les coef beta qui ne sont pas sur la diag et permute
    liste_alpha_permu <- sample(alpha[which(lower.tri(alpha) | upper.tri(alpha))], length(alpha[which(lower.tri(alpha) | upper.tri(alpha))]))## Recup tous les coef beta qui ne sont pas sur la diag et permute
    beta_permu <- beta 
    alpha_permu <- alpha
    c=1
    for (i in 1:S){
      for (j in 1:S){
        if (i!=j){
          beta_permu[i,j] <- liste_beta_permu[c]
          alpha_permu[i,j] <- liste_alpha_permu[c]
          c=c+1
        }
      }
    }
    x = c(x,algo_simul(gen, S, ini, fert,mat_rate,juv_surv,adult_surv,alpha_permu,beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
  }
  return(x)
}