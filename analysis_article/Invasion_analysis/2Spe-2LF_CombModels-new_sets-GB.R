##III - Two-stage model with two-species fertility competition
### Corrected by FB and GB 29/04/2020
### Calculus for all parameter sets of table 2 (Sam boireau's report)
### Parameters are corrected to fit with Ralpha and Rbetas of the table

#1. Clearing workspace
rm(list=ls())
graphics.off()

#2. Parameters set on types (lambda 1, lambda 2, alphas...)

fert1<-c(30,30,30)
fert2<-c(25,25,25)

matur1<-c(0.7,0.7,0.7)
matur2 <-c(0.8,0.8,0.8)

#juvenile survival
phi1 <-c(0.5,0.5,0.5)
phi2 <-c(0.4,0.4,0.4)

#adult survival
s1a<-c(0.5,0.5,0.5)
s2a<-c(0.6,0.6,0.6)

Alphaset<-list(matrix(c(0.1, 0.05, 0.06, 0.1),ncol = 2, byrow = TRUE),
               matrix(c(0.1, 0.02, 0.112, 0.1),ncol = 2, byrow = TRUE),
               matrix(c(0.1, 0.043, 0.035, 0.1),ncol = 2, byrow = TRUE))

Betaset<-list(matrix(c(0.1, 0.06, 0.06, 0.1),ncol = 2, byrow = TRUE),
              matrix(c(0.1, 0.125, 0.01, 0.1),ncol = 2, byrow = TRUE),
              matrix(c(0.1, 0.155, 0.165, 0.1),ncol = 2, byrow = TRUE))


#3. Creating vectors
years <- 3000
TimeVec<-seq(1,years+1)
N1j<-matrix(0,3,years+1)
N2j<-matrix(0,3,years+1)
N1a<-matrix(0,3,years+1)
N2a<-matrix(0,3,years+1)
svar1j<-matrix(0,3,years+1)
svar2j<-matrix(0,3,years+1)

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
    
    svar1j[i,t] <- phi1[i]/(1+betas[1,1]*N1a[i,t]+betas[1,2]*N2a[i,t])
    svar2j[i,t] <- phi2[i]/(1+betas[2,1]*N1a[i,t]+betas[2,2]*N2a[i,t]) ### corrected indices betas[2,1] and betas[2,2] 29/04/2020
    
    N1j[i,t+1]<-(1-matur1[i])*svar1j[i,t]*N1j[i,t]+
      (fert1[i]/(1+alphs[1,1]*N1a[i,t]+alphs[1,2]*N2a[i,t])*N1a[i,t]) ### corrected svar1j[i] -> svar1j[i,t] 29/04/2020
    N2j[i,t+1]<-(1-matur2[i])*svar2j[i,t]*N2j[i,t]+
      (fert2[i]/(1+alphs[2,2]*N2a[i,t]+alphs[2,1]*N1a[i,t])*N2a[i,t]) ### corrected svar2j[i]-> svar2j[i,t] 29/04/2020
    
    N1a[i,t+1]<-matur1[i]*svar1j[i,t]*N1j[i,t]+s1a[i]*N1a[i,t]  ### corrected svar1j[i] -> svar1j[i,t] 29/04/2020
    N2a[i,t+1]<-matur2[i]*svar2j[i,t]*N2j[i,t]+s2a[i]*N2a[i,t]  ### corrected svar2j[i]-> svar2j[i,t] 29/04/2020
  }
}


  #5. Plotting & invasion scores
mat_invasion_crit <-matrix(0,3,4)
par(mfrow=c(2,2))

for(i in 1:3){  

  
  # Invasion scores
  alphs <- Alphaset[[i]]
  betas <- Betaset[[i]]
  
  s1j<- phi1
  s2j <-phi2
  
  scrit1j <- (1-s1a)/((1-matur1)*(1-s1a)+fert1*matur1)
  scrit2j <- (1-s2a)/((1-matur2)*(1-s2a)+fert2*matur2)
  
  f1<-(1/(matur1*s1j))*(1-s1a)*(1-s1j+matur1*s1j)
  f2<-(1/(matur2*s2j))*(1-s2a)*(1-s2j+matur2*s2j)
  
  Rs1<-phi1/scrit1j-1
  Rs2<-phi2/scrit2j-1
  
  Rf1<-(fert1/f1)-1
  Rf2<-(fert2/f2)-1
  
  
  Inv1<-((Rf1[i])/(Rf2[i]))*(alphs[2,2]/alphs[1,2])
  Inv2<-((Rf2[i])/(Rf1[i]))*(alphs[1,1]/alphs[2,1])
  Inv3<-((Rs1[i])/(Rs2[i]))*(betas[2,2]/betas[1,2])
  Inv4<-((Rs2[i])/(Rs1[i]))*(betas[1,1]/betas[2,1])
  
  print(c(Inv1,Inv2,Inv3,Inv4))
  mat_invasion_crit[i,] <- c(Inv1,Inv2,Inv3,Inv4)
  # Plottings
  
  plot(TimeVec,N1j[i,],type="l",
       xlab='Time (years)',ylab='Population Density', cex.main = 1.2, 
       ylim=c(0,1.2*max(N1j[i,], N2j[i,], N1a[i,], N2a[i,])),cex.lab = 1.2,lty = 2, lwd = 2 , col='cyan')
  
  title(paste("Ralpha1=",round(Inv1, 3), " Ralpha2=",round(Inv2, 3), " Rbeta1=", round(Inv3, 3), " Rbeta2=", round(Inv4, 3)),line=0.6,cex=0.7)
  lines(TimeVec,N1a[i,],type="l",lwd=2,col='cyan')
  lines(TimeVec,N2j[i,],lty = 2,lwd=2,col='orange')
  lines(TimeVec,N2a[i,],type="l",lwd=2,col='orange')
  mtext("Combined predictions for the population dynamics of species 1 and species 2", 
        outer=TRUE,  cex=1.2, line=-1.5)
  
  legend("bottomright", legend = c("Species 1", "Species 2"),
         lty = 1, lwd = 2, col = c("cyan", "orange"), cex =1 , bty = "n")
  legend("bottom", legend = c("Juveniles","Adults"),
         lty = c(1,2), cex =1 , bty = "n")
}

N1a[,3000]+N1j[,3000]
N2a[,3000]+N2j[,3000]
