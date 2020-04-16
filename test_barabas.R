n1=0.899
n2=0.114
p=0.36
s=0.3
z=0.7
beta=0.8
F=8
n=n1+n2
R0=2*n2
R1=n1+beta*n1
R2=beta*n1+n1
q=c(1,1)
p_bold=c(n1/n, n2/n)
A1=matrix(nrow=2,ncol=2, rbind(c((1-p)*z/(1+R1), F/(1+R0)),
                     c(p*z/(1+R1), s/(1+R0))))
A2=A1


lambda=eigen(A1)$values
w=eigen(A1)$vectors
w=-1/norm(w, type="O")*w
v=t(eigen(t(A1))$vectors)
v=v*as.numeric(1/(w[,1]%*%v[1,]))

dA1_dR0 = matrix(nrow=2,ncol=2, rbind(c(0, -F/(1+R0)^2),
                                    c(0, -s/(1+R0)^2)))

dA2_dR0 = matrix(nrow=2,ncol=2, rbind(c(0, -F/(1+R0)^2),
                                      c(0, -s/(1+R0)^2)))


dA1_dR1= matrix(nrow=2,ncol=2, rbind(c(-(1-p)*z/(1+R1)^2, 0),
                                 c(-p*z/(1+R1)^2, 0)))

dA1_dR2= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                     c(0, 0)))

dA2_dR1= matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                     c(0, 0)))

dA2_dR2= matrix(nrow=2,ncol=2, rbind(c(-(1-p)*z/(1+R2)^2, 0),
                                     c(-p*z/(1+R2)^2, 0)))

dR0_dn1 = c(0,1)
dR1_dn1=c(1,0)
dR2_dn1=c(beta,0)

dR0_dn2 = c(0,1)
dR1_dn2=c(beta,0)
dR2_dn2=c(1,0)

dG0_dR0 = (n/(q%*%w[,1]))%*%dR0_dn1%*%((1/(lambda[1]-lambda[2]))*kronecker(w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1],v[2,])%*%dA1_dR0%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR0_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR0%*%w[,1])

dG0_dR1 = (n/(q%*%w[,1]))%*%dR0_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR1%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR0_dn2%*%((1/(lambda[1]-lambda[2]))*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR1%*%w[,1])

dG0_dR2 = (n/(q%*%w[,1]))%*%dR0_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR2%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR0_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR2%*%w[,1])

dG1_dR0 = (n/(q%*%w[,1]))%*%dR1_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR0%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR1_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR0%*%w[,1])

dG1_dR1 = (n/(q%*%w[,1]))%*%dR1_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR1%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR1_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR1%*%w[,1])

dG1_dR2 = (n/(q%*%w[,1]))%*%dR1_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR2%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR1_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR2%*%w[,1])

dG2_dR0 = (n/(q%*%w[,1]))%*%dR2_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR0%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR2_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR0%*%w[,1])

dG2_dR1 = (n/(q%*%w[,1]))%*%dR2_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR1%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR2_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR1%*%w[,1])

dG2_dR2 = (n/(q%*%w[,1]))%*%dR2_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dR2%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR2_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dR2%*%w[,1])


dA1_dF = matrix(nrow=2,ncol=2, rbind(c(0, 1/(1+R0)),
                                     c(0, 0)))
dA2_dF = matrix(nrow=2,ncol=2, rbind(c(0, 0),
                                     c(0, 0)))

dG0_dF = (n/(q%*%w[,1]))%*%dR0_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dF%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR0_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dF%*%w[,1])

dG1_dF = (n/(q%*%w[,1]))%*%dR1_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dF%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR1_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dF%*%w[,1])

dG2_dF = (n/(q%*%w[,1]))%*%dR2_dn1%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA1_dF%*%w[,1]) +
  (n/(q%*%w[,1]))%*%dR2_dn2%*%(1/(lambda[1]-lambda[2])*kronecker((w[,2]-((q%*%w[,2])/(q%*%w[,1]))%*%w[,1]),v[2,])%*%dA2_dF%*%w[,1])


dn1_dF = -((v[1,]%*%dA1_dR0%*%w[,1])*((1-dG0_dR0)^-1)*(dR0_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR1)^-1)*(dR1_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR2)^-1)*(dR2_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR0)^-1)*(dR0_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((1-dG1_dR1)^-1)*(dR1_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR2)^-1)*(dR2_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR0)^-1)*(dR0_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR1)^-1)*(dR1_dn1)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((1-dG2_dR2)^-1)*(dR2_dn1)%*%p_bold)^-1 *
  (v[1,]%*%dA1_dF%*%w[,1] + 
    ((v[1,]%*%dA1_dR0%*%w[,1])*((1-dG0_dR0)^-1)* dG0_dF+
     (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR1)^-1)* dG1_dF+
     (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR2)^-1)* dG2_dF+
     (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR0)^-1)* dG0_dF+
     (v[1,]%*%dA1_dR1%*%w[,1])*((1-dG1_dR1)^-1)* dG1_dF+
     (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR2)^-1)* dG2_dF+
     (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR0)^-1)* dG0_dF+
     (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR1)^-1)* dG1_dF+
     (v[1,]%*%dA1_dR2%*%w[,1])*((1-dG2_dR2)^-1)* dG2_dF)) -
  (((v[1,]%*%dA1_dR0%*%w[,1])*((1-dG0_dR0)^-1)*(dR0_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR1)^-1)*(dR1_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR0%*%w[,1])*((0-dG0_dR2)^-1)*(dR2_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR0)^-1)*(dR0_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((1-dG1_dR1)^-1)*(dR1_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR1%*%w[,1])*((0-dG1_dR2)^-1)*(dR2_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR0)^-1)*(dR0_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((0-dG2_dR1)^-1)*(dR1_dn2)%*%p_bold +
  (v[1,]%*%dA1_dR2%*%w[,1])*((1-dG2_dR2)^-1)*(dR2_dn2)%*%p_bold)^-1 *
  (v[1,]%*%dA2_dF%*%w[,1] + 
     ((v[1,]%*%dA2_dR0%*%w[,1])*((1-dG0_dR0)^-1)* dG0_dF+
      (v[1,]%*%dA2_dR0%*%w[,1])*((0-dG0_dR1)^-1)* dG1_dF+
      (v[1,]%*%dA2_dR0%*%w[,1])*((0-dG0_dR2)^-1)* dG2_dF+
      (v[1,]%*%dA2_dR1%*%w[,1])*((0-dG1_dR0)^-1)* dG0_dF+
      (v[1,]%*%dA2_dR1%*%w[,1])*((1-dG1_dR1)^-1)* dG1_dF+
      (v[1,]%*%dA2_dR1%*%w[,1])*((0-dG1_dR2)^-1)* dG2_dF+
      (v[1,]%*%dA2_dR2%*%w[,1])*((0-dG2_dR0)^-1)* dG0_dF+
      (v[1,]%*%dA2_dR2%*%w[,1])*((0-dG2_dR1)^-1)* dG1_dF+
      (v[1,]%*%dA2_dR2%*%w[,1])*((1-dG2_dR2)^-1)* dG2_dF)))



#### TEST EQUILIBRIUM ####

F=8
p=0.36
s=0.3
z=0.7
beta12=0.8
beta21=0.8
n1=c(10,10)
n2=c(10,10)

for (i in 1:100){
  R0=n1[2]+n2[2]
  R1=n1[1]+beta12*n2[1]
  R2=beta21*n1[1]+n2[1]
  A1=matrix(c((1-p)*z/(1+R1), F/(1+R0), p*z/(1+R1), s/(1+R0)), 2,2, byrow = TRUE)
  A2=matrix(c((1-p)*z/(1+R2), F/(1+R0), p*z/(1+R2), s/(1+R0)), 2,2, byrow = TRUE)
  n1=A1%*%n1
  n2=A2%*%n2
}
n1
n2

           