tab_param <- data.frame("pi01"= c(30,30,20,23,30,30,20,30,15,15,10,11,3,30,20,23),
                        "pi02"=c(25,25,25,22,25,25,25,28,12.5,12.5,12.5,11,2.5,25,25,22),
                        "gamma1"=c(0.9,0.7,0.5,0.4,0.9,0.7,0.5,0.4,0.9,0.7,0.5,0.4,0.9,0.7,0.5,0.4),
                        "gamma2"=c(0.8,0.6,0.6,0.4,0.8,0.6,0.6,0.4,0.8,0.6,0.6,0.4,0.8,0.6,0.6,0.4),
                        "phi1"=c(0.5,0.2,0.2,0.19,0.5,0.2,0.2,0.19,0.5,0.2,0.2,0.19,0.5,0.2,0.2,0.19),
                        "phi2"=c(0.4,0.18,0.25,0.2,0.4,0.22,0.25,0.2,0.4,0.18,0.25,0.2,0.4,0.18,0.25,0.2),
                        "s1a"=c(0.5,0.4,0.6,0.2,0.5,0.4,0.4,0.6,0.5,0.4,0.6,0.2,0.5,0.4,0.6,0.2),
                        "s2a"=c(0.6,0.6,0.5,0.4,0.6,0.6,0.3,0.5,0.6,0.6,0.5,0.4,0.6,0.6,0.5,0.4),
                        "alpha11"=c(0.01,0.005,0.01,0.002,0.006,0.01,0.004,0.01,0.01,0.005,0.01,0.002,0.01,0.005,0.01,0.002),
                        "alpha12"=c(0.008,0.004,0.009,0.003,0.007,0.004,0.004,0.005,0.008,0.004,0.009,0.003,0.008,0.004,0.009,0.003),
                        "alpha21"=c(0.0045,0.007,0.0065,0.0045,0.0055,0.006,0.0075,0.0065,0.0045,0.007,0.0065,0.0045,0.0045,0.007,0.0065,0.0045),
                        "alpha22"=c(0.01,0.01,0.01,0.003,0.004,0.01,0.01,0.01,0.01,0.01,0.01,0.003,0.01,0.01,0.01,0.003),
                        "beta11"=c(0.01,0.01,0.003,0.01,0.01,0.004,0.003,0.004,0.01,0.01,0.003,0.01,0.01,0.01,0.003,0.01),
                        "beta12"=c(0.007,0.008,0.003,0.008,0.007,0.003,0.005,0.008,0.007,0.008,0.003,0.008,0.007,0.008,0.003,0.008),
                        "beta21"=c(0.0055,0.004,0.008,0.0045,0.0055,0.007,0.006,0.0065,0.0055,0.004,0.008,0.0045,0.0055,0.004,0.008,0.0045),
                        "beta22"=c(0.01,0.006,0.01,0.004,0.01,0.01,0.004,0.002,0.01,0.006,0.01,0.004,0.01,0.006,0.01,0.004))

tab_param <- as.matrix(tab_param)
result_sensib <- list()
result_inva <- list()

### CALCUL EQUILIBRE ###

# pi01=3
# pi02=2.5
# gamma1 =0.9
# gamma2 =0.8
# phi1=0.5
# phi2=0.4
# s1a=0.5
# s2a =0.6
# n1=c(1,1)
# n2=c(1,1)
# alpha11=0.8
# alpha12=0.7
# alpha21=0.4
# alpha22=0.8
# beta11=0.8
# beta12=0.5
# beta21=0.4
# beta22=0.7
for (k in 1:16){
  k=2
pi01=tab_param[k,1]
pi02=tab_param[k,2]
gamma1=tab_param[k,3]
gamma2=tab_param[k,4]
phi1=tab_param[k,5]
phi2=tab_param[k,6]
s1a=tab_param[k,7]
s2a=tab_param[k,8]
alpha11=tab_param[k,9]
alpha12=tab_param[k,10]
alpha21=tab_param[k,11]
alpha22=tab_param[k,12]
beta11=tab_param[k,13]
beta12=tab_param[k,14]
beta21=tab_param[k,15]
beta22=tab_param[k,16]
n1=c(0,1)
n2=c(0,1)
# ## only fertility ##
# n1=c(1,1)
# n2=c(1,1)
# for (i in 1:1000){
#   R1f=alpha11*n1[2]+alpha12*n2[2]
#   R2f=alpha21*n1[2]+alpha22*n2[2]
#   A1=matrix(c((1-gamma1)*phi1, pi01/(1+R1f), gamma1*phi1, s1a), 2,2, byrow = TRUE)
#   A2=matrix(c((1-gamma2)*phi2, pi02/(1+R2f), gamma2*phi2, s2a), 2,2, byrow = TRUE)
#   n1=A1%*%n1
#   n2=A2%*%n2
# }
# n1
# n2
# Ralpha1 = alpha22*R1f/(R2f*alpha12)
# Ralpha2= alpha11*R2f/(R1f*alpha21)



# ## only survival ##
# n1=c(1,1)
# n2=c(1,1)
# for (i in 1:1000){
#   R1s=beta11*n1[2]+beta12*n2[2]
#   R2s=beta21*n1[2]+beta22*n2[2]
#   A1=matrix(c((1-gamma1)*phi1/(1+R1s), pi01, gamma1*phi1/(1+R1s), s1a), 2,2, byrow = TRUE)
#   A2=matrix(c((1-gamma2)*phi2/(1+R2s), pi02, gamma2*phi2/(1+R2s), s2a), 2,2, byrow = TRUE)
#   n1=A1%*%n1
#   n2=A2%*%n2
# }
# n1
# n2
# Rbeta1 = beta22*R1s/(R2s*beta12)
# Rbeta2= beta11*R2s/(R1s*beta21)


## both ##
n1_both=c(1,1)
n2_both=c(1,1)
for (i in 1:10000){
  R1s=beta11*n1_both[2]+beta12*n2_both[2]
  R2s=beta21*n1_both[2]+beta22*n2_both[2]
  R1f=alpha11*n1_both[2]+alpha12*n2_both[2]
  R2f=alpha21*n1_both[2]+alpha22*n2_both[2]
  A1=matrix(c((1-gamma1)*phi1/(1+R1s), pi01/(1+R1f), gamma1*phi1/(1+R1s), s1a), 2,2, byrow = TRUE)
  A2=matrix(c((1-gamma2)*phi2/(1+R2s), pi02/(1+R2f), gamma2*phi2/(1+R2s), s2a), 2,2, byrow = TRUE)
  n1_both=A1%*%n1_both
  n2_both=A2%*%n2_both
}
n1_both
n2_both
# Ralpha1 = alpha22*R1f/(R2f*alpha12)
# Ralpha2= alpha11*R2f/(R1f*alpha21)
# Rbeta1 = beta22*R1s/(R2f*beta12)
# Rbeta2= beta11*R2f/(R1f*beta21)

R1a<-pi01/((1/(gamma1*phi1))*(1-s1a)*(1-phi1+gamma1*phi1))-1
R2a<-pi02/((1/(gamma2*phi2))*(1-s2a)*(1-phi2+gamma2*phi2))-1
Ralpha1 <- R1a*alpha22/(R2a*alpha12)
Ralpha2 <- R2a*alpha11/(R1a*alpha21)

R1b <- phi1/((1-s1a)/((1-gamma1)*(1-s1a)+pi01*gamma1))-1
R2b <- phi2/((1-s2a)/((1-gamma2)*(1-s2a)+pi02*gamma2))-1

Rbeta1 <- R1b*beta22/(R2b*beta12)
Rbeta2 <- R2b*beta11/(R1b*beta21)
#### ANALYSE SENSIBILITE ####

n1j= n1_both[1]
n1a=n1_both[2]
n2j= n2_both[1]
n2a=n2_both[2]
n1=n1j+n2a
n2=n2j+n2a
q=c(1,1)
p1_bold = c(n1j/n1, n1a/n1)
p2_bold = c(n2j/n1, n2a/n1)

R1s=beta11*n1_both[2]+beta12*n2_both[2]
R2s=beta21*n1_both[2]+beta22*n2_both[2]
R1f=alpha11*n1_both[2]+alpha12*n2_both[2]
R2f=alpha21*n1_both[2]+alpha22*n2_both[2]

A1=matrix(c((1-gamma1)*phi1/(1+R1s), pi01/(1+R1f), gamma1*phi1/(1+R1s), s1a), 2,2, byrow = TRUE)
A2=matrix(c((1-gamma2)*phi2/(1+R2s), pi02/(1+R2f), gamma2*phi2/(1+R2s), s2a), 2,2, byrow = TRUE)

lambda1=eigen(A1)$values
lambda2=eigen(A2)$values

w1=eigen(A1)$vectors
w1=-1/norm(w1, type="O")*w1
v1=eigen(t(A1))$vectors
v1=v1*as.numeric(1/(w1[,1]%*%v1[,1]))

w2=eigen(A1)$vectors
w2=-1/norm(w2, type="O")*w2
v2=eigen(t(A1))$vectors
v2=v2*as.numeric(1/(w2[,1]%*%v2[,1]))

#calcul dA/dR

dA1_dR1s= matrix(nrow=2,ncol=2, rbind(c(-(1-gamma1)*phi1/(1+R1s)^2, 0),
                                      c(-gamma1*phi1/(1+R1s)^2, 0)))
dA1_dR1f= matrix(nrow=2,ncol=2, rbind(c(0, -pi01/(1+R1f)^2),
                                      c(0, 0)))
dA1_dR2s= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                      c(0, 0)))
dA1_dR2f= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                      c(0, 0)))

dA2_dR1s= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                      c(0, 0)))
dA2_dR1f= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                      c(0, 0)))
dA2_dR2s= matrix(nrow=2,ncol=2, rbind(c(-(1-gamma2)*phi2/(1+R2s)^2, 0),
                                      c(-gamma2*phi2/(1+R2s)^2, 0)))
dA2_dR2f= matrix(nrow=2,ncol=2, rbind(c(0, -pi02/(1+R2f)^2),
                                      c(0, 0)))

# calcul des dR/dn
dR1s_dn1 = c(0, beta11)
dR1f_dn1 = c(0, alpha11)
dR2s_dn1 = c(0, beta21)
dR2f_dn1 = c(0, alpha21)

dR1s_dn2 = c(0, beta12) 
dR1f_dn2 = c(0, alpha12)
dR2s_dn2 = c(0, beta22)
dR2f_dn2 = c(0, alpha22)

# calcul dG/dR
prod_kro1 = 1/(lambda1[1]-lambda1[2])*kronecker(matrix(c(w1[,2]-q%*%w1[,2]%*%w1[,1]), 2,1),matrix(v1[,2], 1,2))

prod_kro2 = 1/(lambda2[1]-lambda2[2])*kronecker(matrix(c(w2[,2]-q%*%w2[,2]%*%w2[,1]), 2,1),matrix(v2[,2], 1,2))

dG1s_dR1s= n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
dG1s_dR1f= n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
dG1s_dR2s= n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
dG1s_dR2f= n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]

dG1f_dR1s= n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
dG1f_dR1f= n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
dG1f_dR2s= n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
dG1f_dR2f= n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]

dG2s_dR1s= n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
dG2s_dR1f= n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
dG2s_dR2s= n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
dG2s_dR2f= n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]

dG2f_dR1s= n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
dG2f_dR1f= n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
dG2f_dR2s= n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
dG2f_dR2f= n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]

# calcul dA/dE

dA1_dpi01 = matrix(c(0, 1/(1+R1f), 0, 0), 2,2, byrow=TRUE)
dA2_dpi01 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA1_dpi02 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dpi02 = matrix(c(0, 1/(1+R2f), 0, 0), 2,2, byrow=TRUE)

dA1_dgamma1 = matrix(c(-phi1/(1+R1s), 0, phi1/(1+R1s), 0), 2,2, byrow=TRUE)
dA2_dgamma1 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA1_dgamma2 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dgamma2 = matrix(c(-phi2/(1+R2s), 0, phi2/(1+R2s), 0), 2,2, byrow=TRUE)

dA1_dphi1 = matrix(c((1-gamma1)/(1+R1s), 0, gamma1/(1+R1s), 0), 2,2, byrow=TRUE)
dA2_dphi1 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA1_dphi2 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dphi2 = matrix(c((1-gamma2)/(1+R2s), 0, gamma2/(1+R2s), 0), 2,2, byrow=TRUE)

dA1_ds1a = matrix(c(0, 0, 0, 1), 2,2, byrow=TRUE)
dA2_ds1a = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA1_ds2a = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_ds2a = matrix(c(0, 0, 0, 1), 2,2, byrow=TRUE)

dA1_dalpha11 = matrix(c(0, -n1a*pi01/(1+R1f)^2, 0, 0), 2,2, byrow=TRUE)
dA2_dalpha11 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)

dA1_dalpha12 = matrix(c(0, -n2a*pi01/(1+R1f)^2, 0, 0), 2,2, byrow=TRUE)
dA2_dalpha12 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)

dA1_dalpha21 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dalpha21 = matrix(c(0, -n1a*pi01/(1+R1f)^2, 0, 0), 2,2, byrow=TRUE)

dA1_dalpha22 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dalpha22 = matrix(c(0, -n2a*pi01/(1+R1f)^2, 0, 0), 2,2, byrow=TRUE)

dA1_dbeta11 = matrix(c(-n1a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n1a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)
dA2_dbeta11 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)

dA1_dbeta12 = matrix(c(-n2a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n2a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)
dA2_dbeta12 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)

dA1_dbeta21 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dbeta21 = matrix(c(-n1a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n1a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)

dA1_dbeta22 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
dA2_dbeta22 = matrix(c(-n2a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n1a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)

## calcul dG/dE

dG1s_dpi01 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
dG1s_dpi02 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]

dG1f_dpi01 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
dG1f_dpi02 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]

dG2s_dpi01 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
dG2s_dpi02 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]

dG2f_dpi01 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
dG2f_dpi02 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]

#
dG1s_dgamma1 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
dG1s_dgamma2 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]

dG1f_dgamma1 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
dG1f_dgamma2 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]

dG2s_dgamma1 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
dG2s_dgamma2 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]

dG2f_dgamma1 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
dG2f_dgamma2 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]

#
dG1s_dphi1 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
dG1s_dphi2 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]

dG1f_dphi1 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
dG1f_dphi2 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]

dG2s_dphi1 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
dG2s_dphi2 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]

dG2f_dphi1 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
dG2f_dphi2 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]

#
dG1s_ds1a = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
dG1s_ds2a = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]

dG1f_ds1a = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
dG1f_ds2a = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]

dG2s_ds1a = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
dG2s_ds2a = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]

dG2f_ds1a = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
dG2f_ds2a = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]

#
dG1s_dalpha11 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
dG1s_dalpha12 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
dG1s_dalpha21 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
dG1s_dalpha22 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]

dG1f_dalpha11 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
dG1f_dalpha12 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
dG1f_dalpha21 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
dG1f_dalpha22 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]

dG2s_dalpha11 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
dG2s_dalpha12 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
dG2s_dalpha21 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
dG2s_dalpha22 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]

dG2f_dalpha11 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
dG2f_dalpha12 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
dG2f_dalpha21 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
dG2f_dalpha22 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]

#
dG1s_dbeta11 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
dG1s_dbeta12 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
dG1s_dbeta21 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
dG1s_dbeta22 = n1%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + n2%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]

dG1f_dbeta11 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
dG1f_dbeta12 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
dG1f_dbeta21 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
dG1f_dbeta22 = n1%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + n2%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]

dG2s_dbeta11 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
dG2s_dbeta12 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
dG2s_dbeta21 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
dG2s_dbeta22 = n1%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + n2%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]

dG2f_dbeta11 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
dG2f_dbeta12 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
dG2f_dbeta21 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
dG2f_dbeta22 = n1%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + n2%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]

### Calcul matrice Mij

Mij = matrix(nrow=2,ncol=2,0)
mat_inv = solve((diag(4)-matrix(c(dG1s_dR1s, dG1s_dR1f, dG1s_dR2s, dG1s_dR2f,
                                  dG1f_dR1s, dG1f_dR1f, dG1f_dR2s, dG1f_dR2f,
                                  dG2s_dR1s, dG2s_dR1f, dG2s_dR2s, dG2s_dR2f,
                                  dG2f_dR1s, dG2f_dR1f, dG2f_dR2s, dG2f_dR2f), 4,4, byrow = TRUE))) 

Mij[1,1] = matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%
 mat_inv %*% matrix(c(dR1s_dn1,dR1f_dn1, dR2s_dn1, dR2f_dn1), 4,2, byrow = TRUE)%*%p1_bold

Mij[1,2] = matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%
 mat_inv %*% matrix(c(dR1s_dn2,dR1f_dn2, dR2s_dn2, dR2f_dn2), 4,2, byrow = TRUE)%*%p2_bold

Mij[2,1] = matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%
 mat_inv %*% matrix(c(dR1s_dn1,dR1f_dn1, dR2s_dn1, dR2f_dn1), 4,2, byrow = TRUE)%*%p1_bold

Mij[2,2] = matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%
 mat_inv %*% matrix(c(dR1s_dn2,dR1f_dn2, dR2s_dn2, dR2f_dn2), 4,2, byrow = TRUE)%*%p2_bold

### calcul gj

gj_pi01=c(
  v1[,1]%*%dA1_dpi01%*%w1[,1]+
  matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
  matrix(c(dG1s_dpi01,dG1f_dpi01,dG2s_dpi01,dG2f_dpi01),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dpi01%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dpi01,dG1f_dpi01,dG2s_dpi01,dG2f_dpi01),4,1, byrow = TRUE)
  )
gj_pi02=c(
  v1[,1]%*%dA1_dpi02%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dpi02,dG1f_dpi02,dG2s_dpi02,dG2f_dpi02),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dpi02%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dpi02,dG1f_dpi02,dG2s_dpi02,dG2f_dpi02),4,1, byrow = TRUE)
)
gj_gamma1 =c(
  v1[,1]%*%dA1_dgamma1%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dgamma1,dG1f_dgamma1,dG2s_dgamma1,dG2f_dgamma1),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dgamma1%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dgamma1,dG1f_dgamma1,dG2s_dgamma1,dG2f_dgamma1),4,1, byrow = TRUE)
)
gj_gamma2 =c(
  v1[,1]%*%dA1_dgamma2%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dgamma2,dG1f_dgamma2,dG2s_dgamma2,dG2f_dgamma2),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dgamma2%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dgamma2,dG1f_dgamma2,dG2s_dgamma2,dG2f_dgamma2),4,1, byrow = TRUE)
)
gj_phi1=c(
  v1[,1]%*%dA1_dphi1%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dphi1,dG1f_dphi1,dG2s_dphi1,dG2f_dphi1),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dphi1%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dphi1,dG1f_dphi1,dG2s_dphi1,dG2f_dphi1),4,1, byrow = TRUE)
)
gj_phi2=c(
  v1[,1]%*%dA1_dphi2%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dphi2,dG1f_dphi2,dG2s_dphi2,dG2f_dphi2),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dphi2%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dphi2,dG1f_dphi2,dG2s_dphi2,dG2f_dphi2),4,1, byrow = TRUE)
)
gj_s1a=c(
  v1[,1]%*%dA1_ds1a%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_ds1a,dG1f_ds1a,dG2s_ds1a,dG2f_ds1a),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_ds1a%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_ds1a,dG1f_ds1a,dG2s_ds1a,dG2f_ds1a),4,1, byrow = TRUE)
)
gj_s2a =c(
  v1[,1]%*%dA1_ds2a%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_ds2a,dG1f_ds2a,dG2s_ds2a,dG2f_ds2a),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_ds2a%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_ds2a,dG1f_ds2a,dG2s_ds2a,dG2f_ds2a),4,1, byrow = TRUE)
)
gj_alpha11=c(
  v1[,1]%*%dA1_dalpha11%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha11,dG1f_dalpha11,dG2s_dalpha11,dG2f_dalpha11),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dalpha11%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha11,dG1f_dalpha11,dG2s_dalpha11,dG2f_dalpha11),4,1, byrow = TRUE)
)
gj_alpha12=c(
  v1[,1]%*%dA1_dalpha12%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha12,dG1f_dalpha12,dG2s_dalpha12,dG2f_dalpha12),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dalpha12%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha12,dG1f_dalpha12,dG2s_dalpha12,dG2f_dalpha12),4,1, byrow = TRUE)
)
gj_alpha21=c(
  v1[,1]%*%dA1_dalpha21%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha21,dG1f_dalpha21,dG2s_dalpha21,dG2f_dalpha21),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dalpha21%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha21,dG1f_dalpha21,dG2s_dalpha21,dG2f_dalpha21),4,1, byrow = TRUE)
)
gj_alpha22=c(
  v1[,1]%*%dA1_dalpha22%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha22,dG1f_dalpha22,dG2s_dalpha22,dG2f_dalpha22),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dalpha22%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dalpha22,dG1f_dalpha22,dG2s_dalpha22,dG2f_dalpha22),4,1, byrow = TRUE)
)
gj_beta11=c(
  v1[,1]%*%dA1_dbeta11%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta11,dG1f_dbeta11,dG2s_dbeta11,dG2f_dbeta11),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dbeta11%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta11,dG1f_dbeta11,dG2s_dbeta11,dG2f_dbeta11),4,1, byrow = TRUE)
)
gj_beta12=c(
  v1[,1]%*%dA1_dbeta12%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta12,dG1f_dbeta12,dG2s_dbeta12,dG2f_dbeta12),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dbeta12%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta12,dG1f_dbeta12,dG2s_dbeta12,dG2f_dbeta12),4,1, byrow = TRUE)
)
gj_beta21=c(
  v1[,1]%*%dA1_dbeta21%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta21,dG1f_dbeta21,dG2s_dbeta21,dG2f_dbeta21),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dbeta21%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta21,dG1f_dbeta21,dG2s_dbeta21,dG2f_dbeta21),4,1, byrow = TRUE)
)
gj_beta22=c(
  v1[,1]%*%dA1_dbeta22%*%w1[,1]+
    matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2s%*%w1[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta22,dG1f_dbeta22,dG2s_dbeta22,dG2f_dbeta22),4,1, byrow = TRUE)
  , 
  v2[,1]%*%dA2_dbeta22%*%w2[,1]+
    matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2s%*%w2[,1]), 1,4)%*%mat_inv%*%
    matrix(c(dG1s_dbeta22,dG1f_dbeta22,dG2s_dbeta22,dG2f_dbeta22),4,1, byrow = TRUE)
)

### Calcul sensibilite
Mij_inv = solve(Mij)
dn_dpi01 = -Mij_inv%*%gj_pi01
dn_dpi02 = -Mij_inv%*%gj_pi02
dn_dgamma1 = -Mij_inv%*%gj_gamma1
dn_dgamma2 = -Mij_inv%*%gj_gamma2
dn_dphi1 = -Mij_inv%*%gj_phi1
dn_dphi2 = -Mij_inv%*%gj_phi2
dn_ds1a = -Mij_inv%*%gj_s1a
dn_ds2a = -Mij_inv%*%gj_s2a
dn_dalpha11 = -Mij_inv%*%gj_alpha11
dn_dalpha12 = -Mij_inv%*%gj_alpha12
dn_dalpha21 = -Mij_inv%*%gj_alpha21
dn_dalpha22 = -Mij_inv%*%gj_alpha22
dn_dbeta11 = -Mij_inv%*%gj_beta11
dn_dbeta12 = -Mij_inv%*%gj_beta12
dn_dbeta21 = -Mij_inv%*%gj_beta21
dn_dbeta22 = -Mij_inv%*%gj_beta22
dn1_dE <- c(dn_dpi01[1],dn_dpi02[1],dn_dgamma1[1],dn_dgamma2[1],dn_dphi1[1],dn_dphi2[1],dn_ds1a[1],dn_ds2a[1],dn_dalpha11[1],dn_dalpha12[1],dn_dalpha21[1],dn_dalpha22[1],dn_dbeta11[1],dn_dbeta12[1],dn_dbeta21[1],dn_dbeta22[1])
dn2_dE <- c(dn_dpi01[2],dn_dpi02[2],dn_dgamma1[2],dn_dgamma2[2],dn_dphi1[2],dn_dphi2[2],dn_ds1a[2],dn_ds2a[2],dn_dalpha11[2],dn_dalpha12[2],dn_dalpha21[2],dn_dalpha22[2],dn_dbeta11[2],dn_dbeta12[2],dn_dbeta21[2],dn_dbeta22[2])

result_sensib[[k]]=cbind(tab_param[k,],dn1_dE,dn2_dE)
result_inva[[k]] =c(Ralpha1, Ralpha2, Rbeta1, Rbeta2)
}
View(result_sensib)
View(result_inva)