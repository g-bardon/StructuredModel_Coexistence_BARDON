library("MASS")
library("mvtnorm")


## Fonction pour simuler la dynamique sur 200 generations
## Retourne le nombre d'esp qui subsistent

algo_simul <- function(S, fert, mat_rate, juv_surv, adult_surv, alpha, beta){  
  nj = rep(5, S)
  na = rep(10, S)
  Rs=rep(0, S)
  Rf=rep(0, S)
  A=list()
  for (c in 1:200){
    for (i in 1:S){
      Rs[i] = alpha[,i] %*% na 
      Rf[i] = beta[,i] %*% na 
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)
      n =A[[i]]%*% c(nj[i], na[i])
      nj[i]=n[1]
      na[i]=n[2]
    }
  }
  return(sum(na+nj>1))
}    

## Parametrage avec correlation negative entre aij et bij
parametrage <- function (){
  S = 40
  fert = rlnorm(S, meanlog = 3, sdlog=0.5)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 5,8) ## Loi beta, distribution piquée centrée en 8/13
  juv_surv = rbeta(S, 5,8) ## Loi beta, distribution piquée centrée en 8/13
  adult_surv = rbeta(S, 5,8) ## Loi beta, distribution piquée centrée en 8/13
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.01)
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.003,0.003), sigma=matrix(c(0.001^2, -(0.001^2)*0.9, -(0.001^2)*0.9, 0.001^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  
  c=1
   for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        alpha[i,j] <- alpha_beta_inter[c,1] ## Coef de compet inter dans le reste de la matrice, tiré dans le pool de coef alpha_beta_inter
        beta[i,j] <- alpha_beta_inter[c,2]
        c=c+1
      }
    }
  }
  
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}

param = parametrage()
x_ref = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
beta = param[[7]]


### Permutation des coef pour décorreler les coef alphaij betaij : permuter les coef de beta
x <- c()
for (k in 1:500){  
  liste_beta_permu <- sample(beta[which(lower.tri(beta) | upper.tri(beta))], length(beta[which(lower.tri(beta) | upper.tri(beta))]))  ## Recup tous les coef beta qui ne sont pas sur la diag et permute
  beta_permu <- beta 
  c=1
  for (i in 1:40){
    for (j in 1:40){
      if (i!=j){
        beta_permu[i,j] <- liste_beta_permu[c] ## On remet un coef au hasard sur chaque case de la matrice excepté la diagonale (coef intra)
        c=c+1
      }
    }
  }
  x = c(x,algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}
hist(x, xlim=c(0,40),breaks=1:40)
x_ref
