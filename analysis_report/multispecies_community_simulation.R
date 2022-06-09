library("MASS")
library("mvtnorm")


## Function to simulate the dynamic of the population over 200 generations
## It returns the number of species that subsist at the end of the simulation

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

## Parameterization of to draw the vital rates of species constituting the community -> change parameters of distributions to change the equalization levels
parametrage <- function (){
  S = 40
  fert = rlnorm(S, meanlog = 3, sdlog=0.5)  ## log normale law for the fertility
  mat_rate = rbeta(S, 5,5) ##  beta law
  juv_surv = rbeta(S, 5,5) ##  beta law
  adult_surv = rbeta(S, 5,5) ##  beta law 
  
  alpha_intra = rnorm(S, 0.05, 0.01)  # Normal distribution for the intra-specific coefficients 
  beta_intra = rnorm(S, 0.05, 0.01)
  alpha_beta_inter = rmvnorm(n=S*(S-1), mean=c(0.03,0.03), sigma=matrix(c(0.02^2, -(0.02^2)*0.9, -(0.02^2)*0.9, 0.02^2), 2,2, byrow = TRUE)) ## bivariate normal law for the inter-specific coefficient -> negative correaltion between alpha's and beta's
  alpha=matrix(0,S,S) ## Matrix of all alpha's coefficient :  i colum corresponds to the competition affecting species i
  beta=matrix(0,S,S)  ## Matrix of all beta's coefficient :  i colum corresponds to the competition affecting species i
  
  c=1
   for (i in 1:S){
    for (j in 1:S){
      if (i==j){
        alpha[i,j] <- alpha_intra[i]  ## Intra-specific coeffient on the diagonal
        beta[i,j] <- beta_intra[i]
      }else{
        alpha[i,j] <- alpha_beta_inter[c,1] ## Inter-specific coefficient on the off-diagnonal
        beta[i,j] <- alpha_beta_inter[c,2]
        c=c+1
      }
    }
  }
  
  return(list(S,fert, mat_rate, juv_surv, adult_surv, alpha, beta))   ## Parameter set with correlation between alpha ij and beta ij
}

param = parametrage()
x_ref = algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],param[[7]])
beta = param[[7]]


### Permutation of inter-specific coefficients beta to decorrelate alphaij and betaij
x <- c()
for (k in 1:200){  
  liste_beta_permu <- sample(beta[which(lower.tri(beta) | upper.tri(beta))], length(beta[which(lower.tri(beta) | upper.tri(beta))]))  
  beta_permu <- beta 
  c=1
  for (i in 1:S){
    for (j in 1:S){
      if (i!=j){
        beta_permu[i,j] <- liste_beta_permu[c]
        c=c+1
      }
    }
  }
  x = c(x,algo_simul(param[[1]],param[[2]],param[[3]],param[[4]],param[[5]],param[[6]],beta_permu))  ## index i of x -> result of the dynamic simulation for the parameter set with permutation
}

color= rep("white", S)
color[x_ref] = "cyan"
x<- c(x,x_ref)
hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = "5 species, low equalization")
