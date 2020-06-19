library("MASS")

## sans correlation ##

#plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), 15,10))

simul_sans_correl <- function (){
  S = 20
  fert = rlnorm(S, meanlog = 3, sdlog=0.3)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 15,10) ## Loi beta, distribution piquée centrée en 15/25
  juv_surv = rbeta(S, 15,10) ## Loi beta, distribution piquée centrée en 15/25
  adult_surv = rbeta(S, 15,10) ## Loi beta, distribution piquée centrée en 15/25
  
  alpha_intra = rnorm(S, 0.5, 0.1)  # loi normale. petite variance pour tous les coef de compet, ici parametre ecart-type ( et non variance)
  alpha_inter = rnorm(S, 0.02, 0.003)
  
  beta_intra = rnorm(S, 0.5, 0.1)
  beta_inter = rnorm(S, 0.02, 0.003)
  
  ####
  nj = rep(5, S)
  na = rep(10, S)
  Rs=rep(0, 20)
  Rf=rep(0, 20)
  A=list()
  for (c in 1:100){
    for (i in 1:S){
      Rs[i] = alpha_intra[i] * na[i] + alpha_inter[-i]%*%nj[-i]   ## Calcul des var regulatrice pour chaque espece
      Rf[i] = beta_intra[i] * na[i] + beta_inter[-i]%*%nj[-i]
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)  ## matrice de proj avec les variables regulatrices
      n =A[[i]]%*% c(nj[i], na[i]) ## projection
      nj[i]=n[1] # densite juvenile
      na[i]=n[2] # densite adulte
    }
  }
  return(sum(na+nj>1))   ## Nombre d'espece avec plus d'1 individu
}

## Avec correlation negative entre alphaij et betaij
simul_correl <- function (){
  S = 20
  fert = rlnorm(S, meanlog = 3, sdlog=0.3)
  mat_rate = rbeta(S, 15,10)
  juv_surv = rbeta(S, 15,10)
  adult_surv = rbeta(S, 15,10)
  
  x=rmvnorm(n=20, mean=c(0.02,0.02), sigma=matrix(c(0.003^2,-0.000005,-0.000005,0.003^2), ncol=2))  ### Correlation negative des coef de competion interspécifique, attention ici parametre de la loi = matrice de correl avec les variance/covariance
  alpha_intra = rnorm(S, 0.5, 0.05)
  alpha_inter = x[,1]
  
  beta_intra = rnorm(S, 0.5, 0.05)
  beta_inter = x[,2]
  
  ####
  nj = rep(5, S)
  na = rep(10, S)
  Rs=rep(0, 20)
  Rf=rep(0, 20)
  A=list()
  for (c in 1:100){
    for (i in 1:S){
      Rs[i] = alpha_intra[i] * na[i] + alpha_inter[-i]%*%nj[-i]
      Rf[i] = beta_intra[i] * na[i] + beta_inter[-i]%*%nj[-i]
    }
    for (i in 1:S){
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)   
      n =A[[i]]%*% c(nj[i], na[i])
      nj[i]=n[1]
      na[i]=n[2]
    }
  }
  return(sum(na+nj>1))
}

## Avec correlation negative entre alphaii et betaii

simul_correl2 <- function (){
  S = 20
  fert = rlnorm(S, meanlog = 3, sdlog=0.3)
  mat_rate = rbeta(S, 15,10)
  juv_surv = rbeta(S, 15,10)
  adult_surv = rbeta(S, 15,10)
  
  x=rmvnorm(n=20, mean=c(0.5,0.5), sigma=matrix(c(0.1^2,-0.009,-0.009,0.1^2), ncol=2)) ### Correlation negative des coef de competion intraspécifique
  alpha_intra = x[,1]
  alpha_inter = rnorm(S, 0.02, 0.003)
  
  beta_intra = x[,2]
  beta_inter = rnorm(S, 0.02, 0.003)
  
  ####
  nj = rep(5, S)
  na = rep(10, S)
  Rs=rep(0, 20)
  Rf=rep(0, 20)
  A=list()
  for (c in 1:100){
    for (i in 1:S){
      Rs[i] = alpha_intra[i] * na[i] + alpha_inter[-i]%*%nj[-i]
      Rf[i] = beta_intra[i] * na[i] + beta_inter[-i]%*%nj[-i]
    }
    for (i in 1:S){
      A[[i]]= matrix(c((1-mat_rate[i])*juv_surv[i]/(1+Rs[i]), fert[i]/(1+Rf[i]), mat_rate[i]*juv_surv[i]/(1+Rs[i]), adult_surv[i]), 2,2, byrow = TRUE)
      n =A[[i]]%*% c(nj[i], na[i])
      nj[i]=n[1]
      na[i]=n[2]
    }
  }
  return(sum(na+nj>1))
}


sans_correl=0
avec_correl_inter=0
avec_correl_intra=0
for(i in 1:1000){
  avec_correl_inter = avec_correl_inter + simul_correl()
  avec_correl_intra = avec_correl_intra + simul_correl2()
  sans_correl = sans_correl + simul_sans_correl()
}
avec_correl/1000
sans_correl/1000

## resultat : pas de différence significative entre les différents cas
# -> faire le modèle differement car ici peu d'impact de l'espèce elle même sur sa survie car toute les espèces participent à la competition
# -> Introduire des espèces qui n'interagissent pas 
# -> pour chaque espece, mettre de la correlation positive entre tous les alphaij et entre tous les betaij et une correlation negative entre les alpha et les beta