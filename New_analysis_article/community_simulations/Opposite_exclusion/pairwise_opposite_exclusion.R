### Script to compute the pairwise interaction of species in the case where we permute the competition coefficient initially built to promote pairwise opposite exclusion
## 19/11/2020

library("ggplot2")
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
param <- read.csv("param_vital_rates.csv")
alpha <- as.matrix(read.table("matrix_alpha.csv", sep=" "))
beta <- as.matrix(read.table("matrix_beta.csv", sep=" "))
S=40

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


test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[i,1],param[j,1], param[i,2],param[j,2],param[i,3],param[j,3],param[i,4],param[j,4],alpha_permu[i,i],alpha_permu[j,i],alpha_permu[i,j],alpha_permu[j,j],beta_permu[i,i],beta_permu[j,i],beta_permu[i,j],beta_permu[j,j]))
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

  