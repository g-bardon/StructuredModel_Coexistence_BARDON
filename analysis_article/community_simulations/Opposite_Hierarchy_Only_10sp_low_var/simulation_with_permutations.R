source("simulation_algo.R")

ini <- read.csv("initialization.csv")
param <- read.csv("param_vital_rates.csv")
alpha <- as.matrix(read.table("matrix_alpha.csv", sep=" "))
beta <- as.matrix(read.table("matrix_beta.csv", sep=" "))
gen=3000
S=10

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
  x = c(x,algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha_permu,beta_permu))  ## indice i de x -> result de la dynamique d'un jeu de param avec permutation
}

write.csv(x, file="results_simul_permut.csv", row.names = FALSE)
