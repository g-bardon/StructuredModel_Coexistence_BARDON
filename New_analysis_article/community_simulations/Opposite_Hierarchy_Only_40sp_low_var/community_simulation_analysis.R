source("simulation_algo.R")
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

ini <- read.csv("initialization.csv")
param <- read.csv("param_vital_rates.csv")
alpha <- as.matrix(read.table("matrix_alpha.csv", sep=" "))
beta <- as.matrix(read.table("matrix_beta.csv", sep=" "))
gen=3000
S=40

## Histogram of number of species at equilibrium in original case and permuted cases
x_ref <- algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,beta)
color= rep("white", S)
color[x_ref] = "cyan"
x <- read.csv("results_simul_permut.csv")[,1]
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species, high equalization"))

## 
test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[i,1],param[j,1], param[i,2],param[j,2],param[i,3],param[j,3],param[i,4],param[j,4],alpha[i,i],alpha[j,i],alpha[i,j],alpha[j,j],beta[i,i],beta[j,i],beta[i,j],beta[j,j]))
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

##

cairo_pdf("opposite_hier_dynamic.pdf", width = 10)
traj_ref = algo_simul_traj(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,beta)
plot(1:gen, traj_ref[1,], type="l", ylim = c(0, max(traj_ref)), xlab="Generations", ylab="Individual per species")
for (s in 2:S){
  lines(traj_ref[s,])
}
title("Population dynamics of 40 species with opposite competitive hierarchies")
dev.off()

cairo_ps("opposite_hier_dynamic.eps", width = 10)
traj_ref = algo_simul_traj(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,beta)
plot(1:gen, traj_ref[1,], type="l", ylim = c(0, max(traj_ref)), xlab="Generations", ylab="Individual per species")
for (s in 2:S){
  lines(traj_ref[s,])
}
title("Population dynamics of 40 species with opposite competitive hierarchies")
dev.off()
##

### TEST RANK

liste_beta <- beta[which(lower.tri(beta) | upper.tri(beta))]  
liste_alpha <- alpha[which(lower.tri(alpha) | upper.tri(alpha))]
plot(liste_alpha, liste_beta)
plot(rank(liste_alpha), rank(liste_beta), xlab="Rank of alpha_ij", ylab="Rank of beta_ij")

##

simul_alpha_less <- algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],matrix(0, S, S),beta)
simul_beta_less <- algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,matrix(0, S, S))
