source("initialization_functions.R")

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
for (S in c(5,10,40)){
  print(S)
  for (i in 1:3){
    param = parametrage(S, mu_fert[i], sig_fert[i], ab[i], mu_alpha[i], sig_alpha[i], mu_beta[i], sig_beta[i], sig_inter[i], mu_inter[i])
    ini <- init(S)
    folder = paste(as.character(S), var_list[i], sep="_")
    dir.create(folder)
    df_param_vital_rates <- data.frame("pi"=param[[1]], "gamma"=param[[2]], "phi"=param[[3]], "sa"=param[[4]])
    write.csv(df_param_vital_rates, file = paste(folder,"/param_vital_rates.csv", sep=""), row.names = FALSE)
    write.table(param[[5]], file=paste(folder, "/matrix_alpha.csv", sep=""), row.names = FALSE, col.names = FALSE)
    write.table(param[[6]], file=paste(folder, "/matrix_beta.csv", sep=""), row.names = FALSE, col.names = FALSE)
    write.csv(ini, file = paste(folder,"/initialization.csv", sep=""), row.names = FALSE)
  }
}
