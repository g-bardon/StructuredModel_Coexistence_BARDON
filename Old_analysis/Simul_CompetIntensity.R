library(ggplot2)
library(gridExtra)

### SENSITIVITY ANALYSIS REQUIRES PARAMETERS SET WITH THE FOLLOWING SHAPE : 
#(pi01, pi02, gamma1, gamma2, phi1, phi2, s1a, s2a, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22)

name_param = c("\u03c01", "\u03c02", "\u03b31", "\u03b32", "\u03c61", "\u03c62", "s1a", "s2a", "\u03b111", "\u03b112", "\u03b121", "\u03b122", "\u03b211", "\u03b212", "\u03b221", "\u03b222")

tab=list()
tab_param <- rbind(c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.5,0.1,0.1,0.5,0.5,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.2,0.1,0.1,0.2,0.2,0.1,0.1,0.2),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.11,0.1,0.1,0.11,0.11,0.1,0.1,0.11))


source("sensitivity_analysis_function.R")
ratio = c(5,2,1.1)

n_orig1=c(0,0,0)
n_orig1[1]= colSums(equi_simul(tab_param[1,]))[1]
n_orig1[2]= colSums(equi_simul(tab_param[2,]))[1]
n_orig1[3]= colSums(equi_simul(tab_param[3,]))[1]
n_orig2=c(0,0,0)
n_orig2[1]= colSums(equi_simul(tab_param[1,]))[2]
n_orig2[2]= colSums(equi_simul(tab_param[2,]))[2]
n_orig2[3]= colSums(equi_simul(tab_param[3,]))[2]
n1 = matrix(0,100,3)
n2 = matrix(0,100,3)
for (j in 1:100){
  s_random = runif(1, min=0.5*0.9, max=0.5*1.1)   ### Draw of the new adlut survival in the 10% bandwith around the original value
  tab_param[,7] <- s_random
  n1[j,1] = (colSums(equi_simul(tab_param[1,]))[1] - n_orig1[1])/n_orig1[1]
  n1[j,2] = (colSums(equi_simul(tab_param[2,]))[1] - n_orig1[2])/n_orig1[2]
  n1[j,3] = (colSums(equi_simul(tab_param[3,]))[1] - n_orig1[3])/n_orig1[3]
  
  n2[j,1] = (colSums(equi_simul(tab_param[1,]))[2] - n_orig2[1])/n_orig2[1]
  n2[j,2] = (colSums(equi_simul(tab_param[2,]))[2] - n_orig2[2])/n_orig2[2]
  n2[j,3] = (colSums(equi_simul(tab_param[3,]))[2] - n_orig2[3])/n_orig2[3]
}
n1 <- abs(n1)
n2 <- abs(n2)
colMeans(n1)
colMeans(n2)
