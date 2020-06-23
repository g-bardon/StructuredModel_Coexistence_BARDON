library("MASS")
library("mvtnorm")


## Fonction pour simuler la dynamique sur 100 generations
## Retourne le nombre d'esp qui subsistent

algo_simul <- function(S, fert, mat_rate, juv_surv, adult_surv, alpha, beta){  
nj = rep(5, S)
na = rep(10, S)
Rs=rep(0, 20)
Rf=rep(0, 20)
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
                                                                                   
## Parametrage sans correlation ##

simul_sans_correl <- function (){
  S = 20
  fert = rlnorm(S, meanlog = 2, sdlog=0.3)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 10,10) ## Loi beta, distribution piquée centrée en 0.5
  juv_surv = rbeta(S, 10,10) ## Loi beta, distribution piquée centrée en 0.5
  adult_surv = rbeta(S, 15,10) ## Loi beta, distribution piquée centrée en 15/25
  
  alpha_intra = rnorm(S, 0.5, 0.1)  # Tirage des coef de competition intra selon loi normale 
  alpha_inter = rnorm(S*S, 0.1, 0.005) # Tirage des coef de competition inter selon loi normale : toutes les esp sont en competition avec toutes les autres
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]
      }else{
        alpha[i,j] <- alpha_inter[(i-1)*S + j]
        }
    }
  }
  
  beta_intra = rnorm(S, 0.5, 0.1)
  beta_inter = rnorm(S*S, 0.1, 0.005)
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        beta[i,j] <- beta_intra[i]
      }else{
        beta[i,j] <- beta_inter[(i-1)*S + j]
      }
    }
  }
  
  return(algo_simul(S, fert, mat_rate, juv_surv, adult_surv, alpha, beta))  ## Appel de la fonction de simulation de la dyn 
}

### Parametrage avec chaque espece forte dans l'une des competition et faible dans l'autre 
# Tirage des coefficients selon deux loi normales avec des moyennes differentes

simul_correl <- function(){
  S = 20
  fert = rlnorm(S, meanlog = 2, sdlog=0.3)
  mat_rate = rbeta(S, 10,10)
  juv_surv = rbeta(S, 10,10)
  adult_surv = rbeta(S, 15,10)
  alpha <- matrix(0, 20, 20)
  beta <- matrix(0, 20, 20)
  for (i in 1:S){
    mean_distrib_fort <- rep(0.05, 20)
    mean_distrib_fort[i] <- 0.25
    mean_distrib_faible <- rep(0.15, 20)
    mean_distrib_faible[i] <- 0.75
    
    sigma_distrib_fort <- rep(0.0025^2, 20)
    sigma_distrib_fort[i] <- 0.05^2
    sigma_distrib_faible <- rep(0.0075^2, 20)
    sigma_distrib_faible[i] <- 0.15^2
    
    x1=rmvnorm(n=1, mean=mean_distrib_fort, sigma=diag(sigma_distrib_fort))  ## Tirage des coefficients de compet ou l'esp est favorisee
    x2=rmvnorm(n=1, mean=mean_distrib_faible, sigma=diag(sigma_distrib_faible)) ## Tirage des coefficients de compet ou l'esp est defavorisee
    if (sample(c(0,1), 1)){  ## Tirage si l'esp est favorisee sur les coefficients alpha ou beta (survie juvenile ou fertilité)
      alpha[,i] <- x1
      beta[,i] <- x2
    }else{
      alpha[,i] <- x2
      beta[,i] <- x1
    }
  }
  return(algo_simul(S, fert, mat_rate, juv_surv, adult_surv, alpha, beta))
}



### Comparaison
sans_correl=vector("numeric", 1000)
avec_correl=vector("numeric", 1000)
for(i in 1:1000){
  avec_correl[i] = simul_correl()
  sans_correl[i] = simul_sans_correl()
}
mean(avec_correl)
mean(sans_correl)
hist(avec_correl)
hist(sans_correl)
## resultat : pas de différence significative entre les différents cas
# -> faire le modèle differement car ici peu d'impact de l'espèce elle même sur sa survie car toute les espèces participent à la competition
# -> Introduire des espèces qui n'interagissent pas 
