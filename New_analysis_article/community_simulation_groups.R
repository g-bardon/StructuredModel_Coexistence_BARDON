
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
  for (c in 1:100){
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
  return(na+nj>1)
}    



#### PRIORITY EFFECT + COEXISTENCE
parametrage <- function (S){
  fert = rlnorm(S, meanlog = 3, sdlog=0.005)  ## loi log normale pour la fertilité, assez petite variance sinon très grosse fertilité possibles
  mat_rate = rbeta(S, 500,500) ## Loi beta
  juv_surv = rbeta(S, 500,500) ## Loi beta
  adult_surv = rbeta(S, 500,500) ## Loi beta
  
  alpha_intra = rnorm(S, 0.05, 0.001)  # Tirage des coef de competition intra selon loi normale 
  beta_intra = rnorm(S, 0.05, 0.001)
  sig=0.03
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.04,0.04), sigma=matrix(c(sig^2, -(sig^2)*0.9, -(sig^2)*0.9, sig^2), 2,2, byrow = TRUE)) 
  alpha=matrix(0,S,S)   ## Matrice de tous les coef de competition alpha : colonne i correspond à la competition exercée sur l'espèce i
  beta=matrix(0,S,S)  ## Matrice de tous les coef de competition beta : colonne i correspond à la competition exercée sur l'espèce i
  group=sample(c(0,1), S, replace=TRUE)   ### Each species is either in groupe 0 or group 1
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Coef de compet intra sur la diag
        beta[i,j] <- beta_intra[i]
      }else{
        x <- alpha_beta_inter[c,]
        if (group[i]){ ### If species i is in groupe 1, i is competitive on beta (beta_ij > alpha_ij)
          alpha[i,j] <- min(x) 
          beta[i,j] <- max(x)
          c=c+1
        }else{
          alpha[i,j] <- max(x)
          beta[i,j] <- min(x)
          c=c+1
        }
      }
    }
  }
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Retourne le jeu de parametre avec de la correlation entre alpha ij et beta ij
}


S=40
param = parametrage(S)
alpha= param[[6]]
beta = param[[7]]


x_ref = sum(algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]]))

sp_combined = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
sp_alpha_less = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],matrix(0, 40, 40),param[[7]])
sp_beta_less = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],matrix(0, 40, 40))
df <- data.frame(cbind(sp_alpha_less, sp_beta_less, sp_combined))

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

### TEST RANK

liste_beta <- beta[which(lower.tri(beta) | upper.tri(beta))]  
liste_alpha <- alpha[which(lower.tri(alpha) | upper.tri(alpha))]
plot(liste_alpha, liste_beta)

plot(rank(liste_alpha), rank(liste_beta), xlab="Rank of alpha_ij", ylab="Rank of beta_ij")

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
  x = c(x,sum(algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],alpha_permu,beta_permu)))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}

color= rep("white", S)
color[x_ref] = "cyan"
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species, high equalization"))
