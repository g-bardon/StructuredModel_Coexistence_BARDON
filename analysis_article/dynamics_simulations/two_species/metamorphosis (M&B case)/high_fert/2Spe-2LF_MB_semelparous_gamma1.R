##III - Two-stage model with two-species fertility competition
### Corrected by FB and GB 29/04/2020
### Calculus for all parameter sets of table 2 (Sam boireau's report)
### Parameters are corrected to fit with Ralpha and Rbetas of the table

#1. Clearing workspace
rm(list=ls())
graphics.off()

#2. Parameters set on types (lambda 1, lambda 2, alphas...)

fert1<-c(300,300,300)
fert2<-c(280,280,280)

#matur1<-c(0.7,0.7,0.7)
#matur2 <-c(0.8,0.8,0.8)

matur1<-c(1,1,1)
matur2<-c(1,1,1)

#juvenile survival
phi1 <-c(0.5,0.5,0.5)
phi2 <-c(0.4,0.4,0.4)

#adult survival

#s1a<-c(0.5,0.5,0.5)
#s2a<-c(0.6,0.6,0.6)

s1a<-c(0,0,0)
s2a<-c(0,0,0)

Alphaset<-list(matrix(c(0.1, 0.05, 0.06, 0.1),ncol = 2, byrow = TRUE),
               matrix(c(0.1, 0.02, 0.112, 0.1),ncol = 2, byrow = TRUE),
               matrix(c(0.1, 0.043, 0.035, 0.1),ncol = 2, byrow = TRUE))

Betaset<-list(matrix(c(0.1, 0.06, 0.06, 0.1),ncol = 2, byrow = TRUE),
              matrix(c(0.1, 0.11, 0.05, 0.1),ncol = 2, byrow = TRUE),
              matrix(c(0.1, 0.175, 0.12, 0.1),ncol = 2, byrow = TRUE))


#3. Creating vectors
years <- 3000
TimeVec<-seq(1,years+1)
N1j<-matrix(0,3,years+1)
N2j<-matrix(0,3,years+1)
N1a<-matrix(0,3,years+1)
N2a<-matrix(0,3,years+1)
svar1j<-matrix(0,3,years+1)
svar2j<-matrix(0,3,years+1)

fvar1j<-matrix(0,3,years+1)
fvar2j<-matrix(0,3,years+1)

#4. INITIALIZATION
N0 <- 50
for(i in 1:3){
  N1j[i,1] <-0
  N1a[i,1] <-N0
  N2j[i,1] <-0
  N2a[i,1] <-N0
}

for (t in 1:years){
  for(i in 1:3){ 
    alphs <- Alphaset[[i]]
    betas <- Betaset[[i]]
    
    svar1j[i,t] <- phi1[i]/(1+betas[1,1]*N1j[i,t]+betas[1,2]*N2j[i,t])
    svar2j[i,t] <- phi2[i]/(1+betas[2,1]*N1j[i,t]+betas[2,2]*N2j[i,t]) 
    
    fvar1j[i,t] <- fert1[i]/(1+alphs[1,1]*N1a[i,t]+alphs[1,2]*N2a[i,t])
    fvar2j[i,t] <- fert2[i]/(1+alphs[2,2]*N2a[i,t]+alphs[2,1]*N1a[i,t])
    
    N1j[i,t+1]<-(1-matur1[i])*svar1j[i,t]*N1j[i,t]+ fvar1j[i,t]*N1a[i,t]
    N2j[i,t+1]<-(1-matur2[i])*svar2j[i,t]*N2j[i,t]+ fvar2j[i,t]*N2a[i,t]
    
    N1a[i,t+1]<-matur1[i]*svar1j[i,t]*N1j[i,t]+s1a[i]*N1a[i,t]
    N2a[i,t+1]<-matur2[i]*svar2j[i,t]*N2j[i,t]+s2a[i]*N2a[i,t]
  }
}


pdf(file="2Spe-2LF_MB_semelparous_gamma1.pdf",width = 8, height = 8)
par(mfrow=c(2,2))
for(i in 1:3){  
  plot(TimeVec,N1j[i,]+N1a[i,],type="l",
       xlab='Time (years)',ylab='Population Density', cex.main = 1.2, 
       ylim=c(0,1.2*max(N1j[i,], N2j[i,], N1a[i,], N2a[i,])),cex.lab = 1.2,lty = 2, lwd = 2 , col='cyan')
  lines(TimeVec,N2j[i,]+N2a[i,],lty = 2,lwd=2,col='orange')
}
dev.off()

N1a[,years]+N1j[,years]
N1a[,years-1]+N1j[,years-1]
N2a[,years]+N2j[,years]
N2a[,years-1]+N2j[,years-1]