library("mvtnorm")
init <- function(S){
  nj <- runif(S, min=10, max=50)
  na <- runif(S, min=10, max=50)
  return(cbind(nj,na))
}

parametrage <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.005)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  juv_surv = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  adult_surv = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  
  alpha_intra = rnorm(S, 0.05, 0.001)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.001)
  sig=0.02
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.04,0.04), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        x <- alpha_beta_inter[c,] 
        alpha[i,j] <- x[1] ## Coef de compet inter dans le reste de la matrice, tiré dans le pool de coef alpha_beta_inter
        beta[i,j] <- x[2]
        c=c+1
      }
    }
  }
  
  return(list(fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}

S=10
param = parametrage(S)
ini <- init(S)

df_param_vital_rates <- data.frame("pi"=param[[1]], "gamma"=param[[2]], "phi"=param[[3]], "sa"=param[[4]])
write.csv(df_param_vital_rates, file = "param_vital_rates.csv", row.names = FALSE)
write.table(param[[5]], file="matrix_alpha.csv", row.names = FALSE, col.names = FALSE)
write.table(param[[6]], file="matrix_beta.csv", row.names = FALSE, col.names = FALSE)
write.csv(ini, file = "initialization.csv", row.names = FALSE)
