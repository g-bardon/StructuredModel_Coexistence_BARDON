library("ggplot2")
source("simulation_algo.R")
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
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) #*60/100
  sh <- strheight(text, cex=cex)#*60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

c1<-c("low", "mid", "high")
c2 <- c("5", "10", "40")
c3 <- LETTERS[1:9]
gen=3000
plot_list <- list()
z=1
par(mfrow=c(3,3))
for (i in c1){
  for (j in c2){
    S=as.numeric(j)
    ini <- read.csv(paste("Opposite_Hierarchy_Only_",j,"sp_",i,"_var/initialization.csv", sep = ""))
    param <- read.csv(paste("Opposite_Hierarchy_Only_",j,"sp_",i,"_var/param_vital_rates.csv", sep=""))
    alpha <- as.matrix(read.table(paste("Opposite_Hierarchy_Only_",j,"sp_",i,"_var/matrix_alpha.csv", sep="")))
    beta <- as.matrix(read.table(paste("Opposite_Hierarchy_Only_",j,"sp_",i,"_var/matrix_beta.csv", sep="")))
    
    x_ref <- algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,beta)
    color= rep("white", S)
    color[x_ref] = "gray"
    x <- read.csv(paste("Opposite_Hierarchy_Only_",j,"sp_",i,"_var/results_simul_permut.csv", sep = ""))[,1]
    x<- c(x,x_ref)
    hist(x, xlim=c(0,S),breaks=0:S, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species,",i,"variance"))
    fig_label(c3[z], cex=2)
    z=z+1
  }
}

##
par(mfrow=c(1,1))
alpha <- as.matrix(read.table("Opposite_Hierarchy_Only_40sp_low_var/matrix_alpha.csv"))
beta <- as.matrix(read.table("Opposite_Hierarchy_Only_40sp_low_var/matrix_beta.csv"))
S=as.numeric(j)
ini <- read.csv("Opposite_Hierarchy_Only_40sp_low_var/initialization.csv")
param <- read.csv("Opposite_Hierarchy_Only_40sp_low_var/param_vital_rates.csv")


liste_beta <- beta[which(lower.tri(beta) | upper.tri(beta))]  
liste_alpha <- alpha[which(lower.tri(alpha) | upper.tri(alpha))]
plot(liste_alpha, liste_beta)
plot(rank(liste_alpha), rank(liste_beta), xlab=expression('Rank of' ~ alpha[ij]), ylab=expression('Rank of' ~ beta[ij]))

##
test_inv_score = c()
for (i in 1:S){
  for (j in 1:S){
    if (i != j){
      inv_sc = invasion_score(c(param[i,1],param[j,1], param[i,2],param[j,2],param[i,3],param[j,3],param[i,4],param[j,4],alpha[i,i],alpha[j,i],alpha[i,j],alpha[j,j],beta[i,i],beta[j,i],beta[i,j],beta[j,j]))
      if (inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){ # 1 = coexistence classique
        test_inv_score <- c(test_inv_score, 1)
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]>1 & inv_sc[[3]]<1 & inv_sc[[4]]<1 | inv_sc[[1]]<1 & inv_sc[[2]]<1 & inv_sc[[3]]>1 & inv_sc[[4]]>1){
        test_inv_score <- c(test_inv_score, 2) # priority effect
      }else if(inv_sc[[1]]>1 & inv_sc[[2]]<1 & inv_sc[[3]]<1 & inv_sc[[4]]>1 | inv_sc[[1]]<1 & inv_sc[[2]]>1 & inv_sc[[3]]>1 & inv_sc[[4]]<1){
        test_inv_score <- c(test_inv_score, 3) # coexistence opposite hierarchy
      }else{
        test_inv_score <- c(test_inv_score, 4) ## other cases
      }
    }
  }
}
test_inv_score <- as.factor(test_inv_score)
ggplot(data.frame(test_inv_score), aes(x=test_inv_score)) +
  geom_bar()+
  xlab("1=Coexistence, 2=Priority effect, 3=Opposite exclusion, 4=Other case")+
  labs(title="Outcomes suggested by invasion criteria of pairs of species")