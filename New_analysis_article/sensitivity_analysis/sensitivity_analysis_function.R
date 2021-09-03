### Application of Barabas method (Barabas et al. 2012) to a combined competition model with 2 species and 2 life stages
### Useful function to simulate equilibrium, compute invasion scores and compute sensitivities
### Gael Bardon
### 06/10/2020



equi_simul <- function(parameter){
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
  n1=c(0,30)
  n2=c(0,30)
  lambda_test1=0
  lambda_test2=0
  c=0
  while ((abs(lambda_test1-1)>10^-5 | abs(lambda_test2-1)>10^-5) & c<3000){
    c=c+1
    R1s=beta11*n1[2]+beta12*n2[2]     ### Writing as in Barabas et al. : R1s -> regulating variable associated with competion on juvenile survival
    R2s=beta21*n1[2]+beta22*n2[2]
    R1f=alpha11*n1[2]+alpha12*n2[2]  ###R1f -> regulating variable associated with competition on fertility
    R2f=alpha21*n1[2]+alpha22*n2[2]
    A1=matrix(c((1-gamma1)*phi1/(1+R1s), pi01/(1+R1f), gamma1*phi1/(1+R1s), s1a), 2,2, byrow = TRUE)
    A2=matrix(c((1-gamma2)*phi2/(1+R2s), pi02/(1+R2f), gamma2*phi2/(1+R2s), s2a), 2,2, byrow = TRUE)
    n1=A1%*%n1
    n2=A2%*%n2
    lambda_test1 <- eigen(A1)$value[1]
    lambda_test2 <- eigen(A2)$value[1]
  }
  if (c==3000){
    if (abs(lambda_test1-1)<10^-5){
      return(print("Parameter set does not lead to coexistence and species 1 wins"))
    }else if(abs(lambda_test2-1)<10^-5){
      return(print("Parameter set does not lead to coexistence and species 2 wins"))
    }else{
      return(print("Both species extinct"))
    }
  }
  return(cbind(n1,n2))
} 

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

sensitivity <- function(parameter){
  n <- equi_simul(parameter)
  n1j= n[1,1]
  n1a=n[2,1]
  n2j= n[1,2]
  n2a=n[2,2]
  q=c(1,1)
  n1=as.vector(q%*%c(n1j,n1a))
  n2=as.vector(q%*%c(n2j,n2a))
  p1_bold = c(n1j/n1, n1a/n1)
  p2_bold = c(n2j/n2, n2a/n2)
  
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
  
  R1s=beta11*n1a+beta12*n2a
  R2s=beta21*n1a+beta22*n2a
  R1f=alpha11*n1a+alpha12*n2a
  R2f=alpha21*n1a+alpha22*n2a
  
  A1=matrix(c((1-gamma1)*phi1/(1+R1s), pi01/(1+R1f), gamma1*phi1/(1+R1s), s1a), 2,2, byrow = TRUE)
  A2=matrix(c((1-gamma2)*phi2/(1+R2s), pi02/(1+R2f), gamma2*phi2/(1+R2s), s2a), 2,2, byrow = TRUE)
  
  lambda1=eigen(A1)$values
  lambda2=eigen(A2)$values
  
  w1=eigen(A1)$vectors         ### right eigenvectors
  w1=-1/norm(w1, type="O")*w1
  v1=eigen(t(A1))$vectors  ### left eigenvectors
  v1=v1*as.numeric(1/(w1[,1]%*%v1[,1]))
  
  w2=eigen(A2)$vectors
  w2=-1/norm(w2, type="O")*w2
  v2=eigen(t(A2))$vectors
  v2=v2*as.numeric(1/(w2[,1]%*%v2[,1]))
  
  ### Calculus of partial derivatives of Ai (projection matrix) with respect to regulating variables (R_i,s and R_i,f)
  
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
  
  # Calculus of partial derivatives of regulating variables (R_i,s and R_i,f) with respect to densities vectors n
  
  dR1s_dn1 = c(0, beta11)
  dR1f_dn1 = c(0, alpha11)
  dR2s_dn1 = c(0, beta21)
  dR2f_dn1 = c(0, alpha21)
  
  dR1s_dn2 = c(0, beta12) 
  dR1f_dn2 = c(0, alpha12)
  dR2s_dn2 = c(0, beta22)
  dR2f_dn2 = c(0, alpha22)
  
  # calculus of partial derivatives of G (eq 30 of barabas article) with respect to regulating variables (R_i,s and R_i,f)
  
  prod_kro1 = 1/(lambda1[1]-lambda1[2])*kronecker(matrix(c(w1[,2]-(q%*%w1[,2]/q%*%w1[,1])%*%w1[,1]), 2,1),matrix(v1[,2], 1,2))    ## kronecker product to simplify calculus
  
  prod_kro2 = 1/(lambda2[1]-lambda2[2])*kronecker(matrix(c(w2[,2]-(q%*%w2[,2]/q%*%w2[,1])%*%w2[,1]), 2,1),matrix(v2[,2], 1,2))
  
  dG1s_dR1s= (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
  dG1s_dR1f= (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
  dG1s_dR2s= (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
  dG1s_dR2f= (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]
  
  dG1f_dR1s= (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
  dG1f_dR1f= (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
  dG1f_dR2s= (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
  dG1f_dR2f= (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]
  
  dG2s_dR1s= (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
  dG2s_dR1f= (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
  dG2s_dR2s= (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
  dG2s_dR2f= (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]
  
  dG2f_dR1s= (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dR1s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dR1s%*%w2[,1]
  dG2f_dR1f= (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dR1f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dR1f%*%w2[,1]
  dG2f_dR2s= (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dR2s%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dR2s%*%w2[,1]
  dG2f_dR2f= (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dR2f%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dR2f%*%w2[,1]
  
  ### Calculus of partial derivatives of Ai (projection matrix) with respect to the parameters 
  
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
  dA2_dalpha21 = matrix(c(0, -n1a*pi02/(1+R2f)^2, 0, 0), 2,2, byrow=TRUE)
  
  dA1_dalpha22 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
  dA2_dalpha22 = matrix(c(0, -n2a*pi02/(1+R2f)^2, 0, 0), 2,2, byrow=TRUE)
  
  dA1_dbeta11 = matrix(c(-n1a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n1a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)
  dA2_dbeta11 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
  
  dA1_dbeta12 = matrix(c(-n2a*(1-gamma1)*phi1/(1+R1s)^2, 0, -n2a*gamma1*phi1/(1+R1s)^2, 0), 2,2, byrow=TRUE)
  dA2_dbeta12 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
  
  dA1_dbeta21 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
  dA2_dbeta21 = matrix(c(-n1a*(1-gamma2)*phi2/(1+R2s)^2, 0, -n1a*gamma2*phi2/(1+R2s)^2, 0), 2,2, byrow=TRUE)
  
  dA1_dbeta22 = matrix(c(0, 0, 0, 0), 2,2, byrow=TRUE)
  dA2_dbeta22 = matrix(c(-n2a*(1-gamma2)*phi2/(1+R2s)^2, 0, -n2a*gamma2*phi2/(1+R2s)^2, 0), 2,2, byrow=TRUE)
  
  # calculus of partial derivatives of G (eq 30 of barabas article) with respect to the parameters
  
  dG1s_dpi01 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
  dG1s_dpi02 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]
  
  dG1f_dpi01 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
  dG1f_dpi02 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]
  
  dG2s_dpi01 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
  dG2s_dpi02 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]
  
  dG2f_dpi01 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dpi01%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dpi01%*%w2[,1]
  dG2f_dpi02 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dpi02%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dpi02%*%w2[,1]
  
  #
  dG1s_dgamma1 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
  dG1s_dgamma2 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]
  
  dG1f_dgamma1 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
  dG1f_dgamma2 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]
  
  dG2s_dgamma1 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
  dG2s_dgamma2 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]
  
  dG2f_dgamma1 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dgamma1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dgamma1%*%w2[,1]
  dG2f_dgamma2 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dgamma2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dgamma2%*%w2[,1]
  
  #
  dG1s_dphi1 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
  dG1s_dphi2 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]
  
  dG1f_dphi1 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
  dG1f_dphi2 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]
  
  dG2s_dphi1 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
  dG2s_dphi2 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]
  
  dG2f_dphi1 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dphi1%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dphi1%*%w2[,1]
  dG2f_dphi2 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dphi2%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dphi2%*%w2[,1]
  
  #
  dG1s_ds1a = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
  dG1s_ds2a = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]
  
  dG1f_ds1a = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
  dG1f_ds2a = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]
  
  dG2s_ds1a = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
  dG2s_ds2a = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]
  
  dG2f_ds1a = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_ds1a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_ds1a%*%w2[,1]
  dG2f_ds2a = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_ds2a%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_ds2a%*%w2[,1]
  
  #
  dG1s_dalpha11 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
  dG1s_dalpha12 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
  dG1s_dalpha21 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
  dG1s_dalpha22 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]
  
  dG1f_dalpha11 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
  dG1f_dalpha12 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
  dG1f_dalpha21 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
  dG1f_dalpha22 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]
  
  dG2s_dalpha11 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
  dG2s_dalpha12 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
  dG2s_dalpha21 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
  dG2s_dalpha22 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]
  
  dG2f_dalpha11 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha11%*%w2[,1]
  dG2f_dalpha12 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha12%*%w2[,1]
  dG2f_dalpha21 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha21%*%w2[,1]
  dG2f_dalpha22 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dalpha22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dalpha22%*%w2[,1]
  
  #
  dG1s_dbeta11 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
  dG1s_dbeta12 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
  dG1s_dbeta21 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
  dG1s_dbeta22 = (n1/q%*%w1[,1])%*%dR1s_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1s_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]
  
  dG1f_dbeta11 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
  dG1f_dbeta12 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
  dG1f_dbeta21 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
  dG1f_dbeta22 = (n1/q%*%w1[,1])%*%dR1f_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR1f_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]
  
  dG2s_dbeta11 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
  dG2s_dbeta12 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
  dG2s_dbeta21 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
  dG2s_dbeta22 = (n1/q%*%w1[,1])%*%dR2s_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2s_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]
  
  dG2f_dbeta11 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta11%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta11%*%w2[,1]
  dG2f_dbeta12 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta12%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta12%*%w2[,1]
  dG2f_dbeta21 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta21%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta21%*%w2[,1]
  dG2f_dbeta22 = (n1/q%*%w1[,1])%*%dR2f_dn1%*%prod_kro1%*%dA1_dbeta22%*%w1[,1] + (n2/q%*%w2[,1])%*%dR2f_dn2%*%prod_kro2%*%dA2_dbeta22%*%w2[,1]
  
  ### Calculus of Mij matrix (see eq 32 of barabas et al. 2012)
  
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
  
  ### Calculus of gj vector (see eq 33 of barabas et al. 2012) for each parameter
  
  gj_pi01=c(
    v1[,1]%*%dA1_dpi01%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dpi01,dG1f_dpi01,dG2s_dpi01,dG2f_dpi01),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dpi01%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dpi01,dG1f_dpi01,dG2s_dpi01,dG2f_dpi01),4,1, byrow = TRUE)
  )
  gj_pi02=c(
    v1[,1]%*%dA1_dpi02%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dpi02,dG1f_dpi02,dG2s_dpi02,dG2f_dpi02),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dpi02%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dpi02,dG1f_dpi02,dG2s_dpi02,dG2f_dpi02),4,1, byrow = TRUE)
  )
  gj_gamma1 =c(
    v1[,1]%*%dA1_dgamma1%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dgamma1,dG1f_dgamma1,dG2s_dgamma1,dG2f_dgamma1),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dgamma1%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dgamma1,dG1f_dgamma1,dG2s_dgamma1,dG2f_dgamma1),4,1, byrow = TRUE)
  )
  gj_gamma2 =c(
    v1[,1]%*%dA1_dgamma2%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dgamma2,dG1f_dgamma2,dG2s_dgamma2,dG2f_dgamma2),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dgamma2%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dgamma2,dG1f_dgamma2,dG2s_dgamma2,dG2f_dgamma2),4,1, byrow = TRUE)
  )
  gj_phi1=c(
    v1[,1]%*%dA1_dphi1%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dphi1,dG1f_dphi1,dG2s_dphi1,dG2f_dphi1),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dphi1%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dphi1,dG1f_dphi1,dG2s_dphi1,dG2f_dphi1),4,1, byrow = TRUE)
  )
  gj_phi2=c(
    v1[,1]%*%dA1_dphi2%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dphi2,dG1f_dphi2,dG2s_dphi2,dG2f_dphi2),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dphi2%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dphi2,dG1f_dphi2,dG2s_dphi2,dG2f_dphi2),4,1, byrow = TRUE)
  )
  gj_s1a=c(
    v1[,1]%*%dA1_ds1a%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_ds1a,dG1f_ds1a,dG2s_ds1a,dG2f_ds1a),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_ds1a%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_ds1a,dG1f_ds1a,dG2s_ds1a,dG2f_ds1a),4,1, byrow = TRUE)
  )
  gj_s2a =c(
    v1[,1]%*%dA1_ds2a%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_ds2a,dG1f_ds2a,dG2s_ds2a,dG2f_ds2a),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_ds2a%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_ds2a,dG1f_ds2a,dG2s_ds2a,dG2f_ds2a),4,1, byrow = TRUE)
  )
  gj_alpha11=c(
    v1[,1]%*%dA1_dalpha11%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha11,dG1f_dalpha11,dG2s_dalpha11,dG2f_dalpha11),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dalpha11%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha11,dG1f_dalpha11,dG2s_dalpha11,dG2f_dalpha11),4,1, byrow = TRUE)
  )
  gj_alpha12=c(
    v1[,1]%*%dA1_dalpha12%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha12,dG1f_dalpha12,dG2s_dalpha12,dG2f_dalpha12),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dalpha12%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha12,dG1f_dalpha12,dG2s_dalpha12,dG2f_dalpha12),4,1, byrow = TRUE)
  )
  gj_alpha21=c(
    v1[,1]%*%dA1_dalpha21%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha21,dG1f_dalpha21,dG2s_dalpha21,dG2f_dalpha21),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dalpha21%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha21,dG1f_dalpha21,dG2s_dalpha21,dG2f_dalpha21),4,1, byrow = TRUE)
  )
  gj_alpha22=c(
    v1[,1]%*%dA1_dalpha22%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha22,dG1f_dalpha22,dG2s_dalpha22,dG2f_dalpha22),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dalpha22%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dalpha22,dG1f_dalpha22,dG2s_dalpha22,dG2f_dalpha22),4,1, byrow = TRUE)
  )
  gj_beta11=c(
    v1[,1]%*%dA1_dbeta11%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta11,dG1f_dbeta11,dG2s_dbeta11,dG2f_dbeta11),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dbeta11%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta11,dG1f_dbeta11,dG2s_dbeta11,dG2f_dbeta11),4,1, byrow = TRUE)
  )
  gj_beta12=c(
    v1[,1]%*%dA1_dbeta12%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta12,dG1f_dbeta12,dG2s_dbeta12,dG2f_dbeta12),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dbeta12%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta12,dG1f_dbeta12,dG2s_dbeta12,dG2f_dbeta12),4,1, byrow = TRUE)
  )
  gj_beta21=c(
    v1[,1]%*%dA1_dbeta21%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta21,dG1f_dbeta21,dG2s_dbeta21,dG2f_dbeta21),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dbeta21%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta21,dG1f_dbeta21,dG2s_dbeta21,dG2f_dbeta21),4,1, byrow = TRUE)
  )
  gj_beta22=c(
    v1[,1]%*%dA1_dbeta22%*%w1[,1]+
      matrix(c(v1[,1]%*%dA1_dR1s%*%w1[,1],v1[,1]%*%dA1_dR1f%*%w1[,1],v1[,1]%*%dA1_dR2s%*%w1[,1], v1[,1]%*%dA1_dR2f%*%w1[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta22,dG1f_dbeta22,dG2s_dbeta22,dG2f_dbeta22),4,1, byrow = TRUE)
    , 
    v2[,1]%*%dA2_dbeta22%*%w2[,1]+
      matrix(c(v2[,1]%*%dA2_dR1s%*%w2[,1],v2[,1]%*%dA2_dR1f%*%w2[,1],v2[,1]%*%dA2_dR2s%*%w2[,1], v2[,1]%*%dA2_dR2f%*%w2[,1]), 1,4)%*%mat_inv%*%
      matrix(c(dG1s_dbeta22,dG1f_dbeta22,dG2s_dbeta22,dG2f_dbeta22),4,1, byrow = TRUE)
  )
  
  ### Calculus of sensitivities (see eq 31)
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
  elasticity_table <- cbind(parameter, (parameter/n1)*dn1_dE,(parameter/n2)*dn2_dE)
  sensitivity_table <- cbind(parameter, dn1_dE, dn2_dE)
  return(list("Mij_matrix"=Mij, "Determinant"=det(Mij), "Sensitivity"=sensitivity_table, "Elasticity"=elasticity_table))
}

### TEST ###
equi_simul(c(30,25,  
             0.8,0.8,
             0.5,0.4,
             0.4,0.6,
             0.1,0.08,0.08,0.1,
             0.1,0.08,0.08,0.1))
invasion_score(c(30,25,  
                 0.8,0.8,
                 0.5,0.4,
                 0.4,0.6,
                 0.1,0.08,0.08,0.1,
                 0.1,0.08,0.08,0.1))

sensitivity(c(30,25,  
              0.8,0.8,
              0.5,0.4,
              0.4,0.6,
              0.1,0.08,0.08,0.1,
              0.1,0.08,0.08,0.1))
