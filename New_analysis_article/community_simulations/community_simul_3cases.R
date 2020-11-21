#Author : Gaël Bardon
#Date : 06/10/2020

# Script to simulate dynamics of community in three situations :
# - for each pair of species, priority effect suggested by beta and coexistence suggested by alpha
# - for each pair of species, priority effect and coexistence suggested OR opposite exclusion
# - for each pair of species, only opposite exclusion

library("MASS")
library("mvtnorm")
library(ggplot2)

## Calcul critère invasion paire
invasion_score <- function(parameter){
  pi01=parameter[1]
  pi02=parameter[2]
  gamma1=parameter[3]
  gamma2=parameter[4]
  phi1=parameter[5]
  phi2=parameter[6]
  s1a=parameter[7]
  s2a=parameter[8]
  alpha11=parameter[9]
  alpha12=parameter[10]
  alpha21=parameter[11]
  alpha22=parameter[12]
  beta11=parameter[13]
  beta12=parameter[14]
  beta21=parameter[15]
  beta22=parameter[16]
  
  scrit1j <- (1-s1a)/((1-gamma1)*(1-s1a)+pi01*gamma1)
  scrit2j <- (1-s2a)/((1-gamma2)*(1-s2a)+pi02*gamma2)
  
  f1<-(1/(gamma1*phi1))*(1-s1a)*(1-phi1+gamma1*phi1)
  f2<-(1/(gamma2*phi2))*(1-s2a)*(1-phi2+gamma2*phi2)
  
  Rs1<-phi1/scrit1j-1
  Rs2<-phi2/scrit2j-1
  
  Rf1<-(pi01/f1)-1
  Rf2<-(pi02/f2)-1
  
  
  Inv1<-(Rf1/Rf2)*(alpha22/alpha12)
  Inv2<-(Rf2/Rf1)*(alpha11/alpha21)
  Inv3<-(Rs1/Rs2)*(beta22/beta12)
  Inv4<-(Rs2/Rs1)*(beta11/beta21)
  
  return(list("Ralpha1"=Inv1,"Ralpha2"=Inv2,"Rbeta1"=Inv3,"Rbeta2"=Inv4))
}


## Fonction pour simuler la dynamique sur 200 generations
## Retourne le nombre d'esp qui subsistent

algo_simul <- function(S, fert, mat_rate, juv_surv, adult_surv, alpha, beta){  
  nj = rep(5, S)
  na = rep(10, S)
  nj_plus = rep(0, S)
  na_plus = rep(0, S)
  Rs=rep(0, S)
  Rf=rep(0, S)
  A=list()
  for (c in 1:3000){
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

#### PRIORITY EFFECT
parametrage_PE <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.05)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  juv_surv = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  adult_surv = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.01)
  sig=0.005
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.01,0.09), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        x <- alpha_beta_inter[c,] ## only priority effect
        alpha[i,j] <- x[1] ## Coef de compet inter dans le reste de la matrice, tiré dans le pool de coef alpha_beta_inter
        beta[i,j] <- x[2]
        c=c+1
      }
    }
  }
  
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}

S=40
param = parametrage_PE(S)

x_ref = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
alpha= param[[6]]
beta = param[[7]]

test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[[2]][i],param[[2]][j], param[[3]][i],param[[3]][j],param[[4]][i],param[[4]][j],param[[5]][i],param[[5]][j],param[[6]][i,i],param[[6]][j,i],param[[6]][i,j],param[[6]][j,j],param[[7]][i,i],param[[7]][j,i],param[[7]][i,j],param[[7]][j,j]))
      if (inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){ # 1 = coexistence classique
        test_inv_score <- c(test_inv_score, 1)
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]<1 & inv_sc[[4]]<1 | inv_sc[[1]]<1 & inv_sc[[2]]<1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){
        test_inv_score <- c(test_inv_score, 2) # priority effect
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]<1 & inv_sc[[3]]<1 & inv_sc[[4]]>1 | inv_sc[[1]]<1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]<1){
        test_inv_score <- c(test_inv_score, 3) # coexistence opposite hierarchy
      }else{
        test_inv_score <- c(test_inv_score, 4) ## other cases
      }
    }
  }
}
test_inv_score <- as.factor(test_inv_score)

ggplot(data.frame(test_inv_score), aes(x=test_inv_score)) +
  geom_bar()+
  xlab("1=Coexistence, 2=Priority effect, 3=Opposite exclusion, 4=Other case")+
  labs(title="Outcomes suggested by invasion criteria of pairs of species")

### Permutation within all inter-specific competition coefficient alpha and within all inter-specific competition coefficient beta

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
        alpha_permu[i,j] <- liste_alpha_permu[c]## On remet un coef au hasard sur chaque case de la matrice excepté la diagonale (coef intra)
        c=c+1
      }
    }
  }
  x = c(x,algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],alpha_permu,beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}

color= rep("white", S)
color[x_ref] = "cyan"
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species, high equalization"))

### PRIORITY EFFECT + COEXISTENCE OR OPPOSITE EXCLUSION

parametrage_PE_OE <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.05)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  juv_surv = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  adult_surv = rbeta(S, 50,50) ## Loi beta, distribution piquée centrée en 8/13
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.01)
  sig=0.005
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.01,0.09), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        x <- sample(alpha_beta_inter[c,]) ## priority effect + opposite exclusion
        alpha[i,j] <- x[1] ## Coef de compet inter dans le reste de la matrice, tiré dans le pool de coef alpha_beta_inter
        beta[i,j] <- x[2]
        c=c+1
      }
    }
  }
  
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}

S=40
param = parametrage_PE_OE(S)
#param[[6]][which(param[[6]]< 0)] <- 0.001
#param[[7]][which(param[[7]]< 0)] <- 0.001
x_ref = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
alpha= param[[6]]
beta = param[[7]]

test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[[2]][i],param[[2]][j], param[[3]][i],param[[3]][j],param[[4]][i],param[[4]][j],param[[5]][i],param[[5]][j],param[[6]][i,i],param[[6]][j,i],param[[6]][i,j],param[[6]][j,j],param[[7]][i,i],param[[7]][j,i],param[[7]][i,j],param[[7]][j,j]))
      if (inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){ # 1 = coexistence classique
        test_inv_score <- c(test_inv_score, 1)
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]<1 & inv_sc[[4]]<1 | inv_sc[[1]]<1 & inv_sc[[2]]<1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){
        test_inv_score <- c(test_inv_score, 2) # priority effect
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]<1 & inv_sc[[3]]<1 & inv_sc[[4]]>1 | inv_sc[[1]]<1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]<1){
        test_inv_score <- c(test_inv_score, 3) # coexistence opposite hierarchy
      }else{
        test_inv_score <- c(test_inv_score, 4) ## other cases
      }
    }
  }
}
test_inv_score <- as.factor(test_inv_score)

ggplot(data.frame(test_inv_score), aes(x=test_inv_score)) +
  geom_bar()+
  xlab("1=Coexistence, 2=Priority effect, 3=Opposite exclusion, 4=Other case")+
  labs(title="Outcomes suggested by invasion criteria of pairs of species")

### Permutation des coef pour décorreler les coef alphaij betaij : permuter les coef de beta

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
        alpha_permu[i,j] <- liste_alpha_permu[c]## On remet un coef au hasard sur chaque case de la matrice excepté la diagonale (coef intra)
        c=c+1
      }
    }
  }
  x = c(x,algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],alpha_permu,beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}

color= rep("white", S)
color[x_ref] = "cyan"
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species, high equalization"))

### OPPOSITE EXCLUSION

parametrage_OE <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.005)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  juv_surv = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  adult_surv = rbeta(S, 500,500) ## Loi beta, distribution piquée centrée en 8/13
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.01)
  sig=0.005
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.01,0.09), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  
  c=1
  for (i in 1:S){
    for (j in 1:i){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        x1 <- alpha_beta_inter[c,] ## priority effect + opposite exclusion
        x2 <- alpha_beta_inter[c+1,]
        alpha[i,j] <- max(x1) ## Coef de compet inter dans le reste de la matrice, tiré dans le pool de coef alpha_beta_inter
        beta[i,j] <- min(x1)
        alpha[j,i] <- min(x2)
        beta[j,i] <- max(x2)
        c=c+2
      }
    }
  }
  
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}

S=40
param = parametrage_OE(S)

x_ref = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
alpha= param[[6]]
beta = param[[7]]

liste_beta <- beta[which(lower.tri(beta) | upper.tri(beta))]  
liste_alpha <- alpha[which(lower.tri(alpha) | upper.tri(alpha))]
plot(liste_alpha, liste_beta)
plot(rank(liste_alpha), rank(liste_beta), xlab="Rank of alpha_ij", ylab="Rank of beta_ij")

test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[[2]][i],param[[2]][j], param[[3]][i],param[[3]][j],param[[4]][i],param[[4]][j],param[[5]][i],param[[5]][j],param[[6]][i,i],param[[6]][j,i],param[[6]][i,j],param[[6]][j,j],param[[7]][i,i],param[[7]][j,i],param[[7]][i,j],param[[7]][j,j]))
      if (inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){ # 1 = coexistence classique
        test_inv_score <- c(test_inv_score, 1)
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]<1 & inv_sc[[4]]<1 | inv_sc[[1]]<1 & inv_sc[[2]]<1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){
        test_inv_score <- c(test_inv_score, 2) # priority effect
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]<1 & inv_sc[[3]]<1 & inv_sc[[4]]>1 | inv_sc[[1]]<1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]<1){
        test_inv_score <- c(test_inv_score, 3) # coexistence opposite hierarchy
      }else{
        test_inv_score <- c(test_inv_score, 4) ## other cases
      }
    }
  }
}
test_inv_score <- as.factor(test_inv_score)

ggplot(data.frame(test_inv_score), aes(x=test_inv_score)) +
  geom_bar()+
  xlab("1=Coexistence, 2=Priority effect, 3=Opposite exclusion, 4=Other case")+
  labs(title="Outcomes suggested by invasion criteria of pairs of species")

### Permutation des coef pour décorreler les coef alphaij betaij : permuter les coef de beta

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
        alpha_permu[i,j] <- liste_alpha_permu[c]  ## On remet un coef au hasard sur chaque case de la matrice excepté la diagonale (coef intra)
        c=c+1
      }
    }
  }
  x = c(x,algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],alpha_permu,beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}

color= rep("white", S)
color[x_ref] = "gray"
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species, low variance"))

### Invasion test on permuted parameter set
test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[[2]][i],param[[2]][j], param[[3]][i],param[[3]][j],param[[4]][i],param[[4]][j],param[[5]][i],param[[5]][j],alpha_permu[i,i],alpha_permu[j,i],alpha_permu[i,j],alpha_permu[j,j],beta_permu[i,i],beta_permu[j,i],beta_permu[i,j],beta_permu[j,j]))
      if (inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){ # 1 = coexistence classique
        test_inv_score <- c(test_inv_score, 1)
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]<1 & inv_sc[[4]]<1 | inv_sc[[1]]<1 & inv_sc[[2]]<1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){
        test_inv_score <- c(test_inv_score, 2) # priority effect
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]<1 & inv_sc[[3]]<1 & inv_sc[[4]]>1 | inv_sc[[1]]<1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]<1){
        test_inv_score <- c(test_inv_score, 3) # coexistence opposite hierarchy
      }else{
        test_inv_score <- c(test_inv_score, 4) ## other cases
      }
    }
  }
}
test_inv_score <- as.factor(test_inv_score)
ggplot(data.frame(test_inv_score), aes(x=test_inv_score)) +
  geom_bar()+
  xlab("1=Coexistence, 2=Priority effect, 3=Opposite exclusion, 4=Other case")+
  labs(title="Outcomes suggested by invasion criteria of pairs of species")

liste_beta_permu <- beta_permu[which(lower.tri(beta_permu) | upper.tri(beta_permu))]  
liste_alpha_permu <- alpha_permu[which(lower.tri(alpha_permu) | upper.tri(alpha_permu))]
plot(liste_alpha_permu, liste_beta_permu)
plot(rank(liste_alpha_permu), rank(liste_beta_permu), xlab="Rank of alpha_ij", ylab="Rank of beta_ij")
