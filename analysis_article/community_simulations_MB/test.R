mu_fert = c(3,3,3)
sig_fert=c(0.005,0.05,0.5)
ab= c(500,50,5)
mu_alpha=c(0.05,0.05,0.05)
sig_alpha=c(0.001,0.01,0.01)
mu_beta=c(0.05,0.05,0.05)
sig_beta=c(0.001,0.01,0.01)
sig_inter=c(0.02,0.02,0.02)
mu_inter=c(0.04,0.04,0.04)
var_list = c("low", "average", "high")
source("initialization_functions.R")
source("simul_functions.R")

S=40
x=c()
x_ref=c()
gen=3000
for (i in 1:300){
  param = parametrage(S, mu_fert[1], sig_fert[1], ab[1], mu_alpha[1], sig_alpha[1], mu_beta[1], sig_beta[1], sig_inter[1], mu_inter[1])
  ini <- init(S)
  x_ref <- c(x_ref, algo_simul(gen, S, ini, param[[1]] ,param[[2]],param[[3]],param[[4]],param[[5]],param[[6]]))
  x <-c(x,simul_permut(gen, S, ini, param[[1]] ,param[[2]],param[[3]],param[[4]], param[[5]],param[[6]]))
  
}

hist(x)
hist(x_ref)