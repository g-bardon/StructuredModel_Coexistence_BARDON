library(ggplot2)
library(gridExtra)
library(Cairo)
### SENSITIVITY ANALYSIS REQUIRES PARAMETERS SET WITH THE FOLLOWING SHAPE : 
#(pi01, pi02, gamma1, gamma2, phi1, phi2, s1a, s2a, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22)
name_param = c("pi1", "pi2", "gamma1", "gamma2", "phi1", "phi2", "s1a", "s2a", "alpha11", "alpha12", "alpha21", "alpha22", "beta11", "beta12", "beta21", "beta22")

tab_param <- rbind(c(30,25,0.7,0.8,0.5,0.4,0.5,0.6, 0.1, 0.05, 0.06, 0.1,0.1, 0.06, 0.06, 0.1),
                    c(30,25,0.7,0.8,0.5,0.4,0.5,0.6, 0.1, 0.02, 0.112, 0.1,0.1, 0.125, 0.01, 0.1),
                   c(30,25,0.7,0.8,0.5,0.4,0.5,0.6, 0.1, 0.043, 0.035, 0.1,0.1, 0.155, 0.165, 0.1))

source("C:/Users/Gael/Documents/Stage_IMB/Scripts R/New_analysis_article/sensitivity_analysis_function.R")
equi_simul(tab_param[1,])
invasion_score(tab_param[1,])
sensitivity(tab_param[1,])

equi_simul(tab_param[2,])
invasion_score(tab_param[2,])
sensitivity(tab_param[2,])

equi_simul(tab_param[3,])
invasion_score(tab_param[3,])
sensitivity(tab_param[3,])



elas1 <- sensitivity(tab_param[1,])$Elasticity
elas2 <- sensitivity(tab_param[2,])$Elasticity
elas3 <- sensitivity(tab_param[3,])$Elasticity

sensi1 <- sensitivity(tab_param[1,])$Sensitivity
sensi2 <- sensitivity(tab_param[2,])$Sensitivity
sensi3 <- sensitivity(tab_param[3,])$Sensitivity

####

sensi1 <- cbind(sensi1, 1)
sensi2 <- cbind(sensi2, 2)
sensi3 <- cbind(sensi3, 3)
data <- rbind(sensi1, sensi2, sensi3)
data[,1] <- c(name_param, name_param, name_param)
data <- as.data.frame(data)
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))

# Grouped species 1
ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + scale_fill_grey() +
  scale_x_discrete(limits = name_param, labels= 
                     c("pi1"=expression(pi[1]), "pi2"=expression(pi[2]), "gamma1"=expression(gamma[1]), "gamma2"=expression(gamma[2]), 
                       "phi1"=expression(phi[1]), "phi2"=expression(phi[2]), "s1a"=expression("s"["1a"]), "s2a"=expression("s"["2a"]),
                       "alpha11"=expression(alpha[11]), "alpha12"=expression(alpha[12]), "alpha21"=expression(alpha[21]), "alpha22"=expression(alpha[22]),
                       "beta11"=expression(beta[11]), "beta12"=expression(beta[12]), "beta21"=expression(beta[21]), "beta22"=expression(beta[22]))) +
  ylab("Sensitivity") + xlab("") +
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Sensitivities of the parameters for the species 1")


# Grouped species 2
ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + scale_fill_grey() +
  scale_x_discrete(limits = name_param, labels= 
                     c("pi1"=expression(pi[1]), "pi2"=expression(pi[2]), "gamma1"=expression(gamma[1]), "gamma2"=expression(gamma[2]), 
                       "phi1"=expression(phi[1]), "phi2"=expression(phi[2]), "s1a"=expression("s"["1a"]), "s2a"=expression("s"["2a"]),
                       "alpha11"=expression(alpha[11]), "alpha12"=expression(alpha[12]), "alpha21"=expression(alpha[21]), "alpha22"=expression(alpha[22]),
                       "beta11"=expression(beta[11]), "beta12"=expression(beta[12]), "beta21"=expression(beta[21]), "beta22"=expression(beta[22]))) +
  ylab("Sensitivity") + xlab("") + 
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Sensitivities of the parameters for the species 2")
s<-expression(pi)
####

elas1 <- cbind(elas1, 1)
elas2 <- cbind(elas2, 2)
elas3 <- cbind(elas3, 3)
data <- rbind(elas1, elas2, elas3)
data[,1] <- c(name_param, name_param, name_param)
data <- as.data.frame(data)
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))
# Grouped species 1
elast_sp1 = ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  +  scale_fill_grey() +
  scale_x_discrete(limits = name_param, labels= 
                     c("pi1"=expression(pi[1]), "pi2"=expression(pi[2]), "gamma1"=expression(gamma[1]), "gamma2"=expression(gamma[2]), 
                       "phi1"=expression(phi[1]), "phi2"=expression(phi[2]), "s1a"=expression("s"["1a"]), "s2a"=expression("s"["2a"]),
                       "alpha11"=expression(alpha[11]), "alpha12"=expression(alpha[12]), "alpha21"=expression(alpha[21]), "alpha22"=expression(alpha[22]),
                       "beta11"=expression(beta[11]), "beta12"=expression(beta[12]), "beta21"=expression(beta[21]), "beta22"=expression(beta[22]))) +
  ylab("Elasticity") + xlab("") +
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Elasticities of the parameters for the species 1")

# Grouped species 2
elast_sp2 = ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + scale_fill_grey() +
  scale_x_discrete(limits = name_param, labels= 
                     c("pi1"=expression(pi[1]), "pi2"=expression(pi[2]), "gamma1"=expression(gamma[1]), "gamma2"=expression(gamma[2]), 
                       "phi1"=expression(phi[1]), "phi2"=expression(phi[2]), "s1a"=expression("s"["1a"]), "s2a"=expression("s"["2a"]),
                       "alpha11"=expression(alpha[11]), "alpha12"=expression(alpha[12]), "alpha21"=expression(alpha[21]), "alpha22"=expression(alpha[22]),
                       "beta11"=expression(beta[11]), "beta12"=expression(beta[12]), "beta21"=expression(beta[21]), "beta22"=expression(beta[22]))) +
  
  ylab("Elasticity") + xlab("") + 
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Elasticities of the parameters for the species 2")

cairo_pdf("combined_elast_3sets.pdf")
plot_grid(elast_sp1, elast_sp2, ncol=1, nrow=2)
dev.off()

cairo_ps("combined_elast_3sets.eps")
plot_grid(elast_sp1, elast_sp2, ncol=1, nrow=2)
dev.off()

a=list()
param_set_index = c(1,2,3)
for (i in 1:3){
  n=colSums(equi_simul(tab_param[i,]))
  coex_region <- sensitivity(tab_param[i,])$Sensitivity
  borne_coexist <- cbind(coex_region[,1], -n[1]/coex_region[,2] , -n[2]/coex_region[,3] )
  for (j in 1:nrow(borne_coexist)){
    borne_coexist[j,2:3] <- c(min(borne_coexist[j,2:3]), max(borne_coexist[j,2:3]))
  }
  borne_coexist = as.data.frame(borne_coexist)
  rownames(borne_coexist) <- c("pi1", "pi2", "gamma1", "gamma2", "phi1", "phi2", "s1a", "s2a", "alpha11", "alpha12", "alpha21", "alpha22", "beta11", "beta12", "beta21", "beta22")
  plot1 <- ggplot(borne_coexist[1:2,], aes(x=rownames(borne_coexist[1:2,]), y=borne_coexist[1:2,1])) + 
    geom_pointrange(ymin=(borne_coexist[1:2,1]+borne_coexist[1:2,2]), ymax=(borne_coexist[1:2,1]+borne_coexist[1:2,3]))+ 
    coord_cartesian(ylim = c(0, 50))+xlab("Fertility") + ylab("Parameter value")+
    scale_x_discrete(labels= c("pi1"=expression(pi[1]), "pi2"=expression(pi[2])))
  plot2 <- ggplot(borne_coexist[3:8,], aes(x=rownames(borne_coexist[3:8,]), y=borne_coexist[3:8,1])) + 
    geom_pointrange(ymin=(borne_coexist[3:8,1]+borne_coexist[3:8,2]), ymax=(borne_coexist[3:8,1]+borne_coexist[3:8,3]))+ 
    coord_cartesian(ylim = c(0, 1))+
    xlab("Maturation rate, juvenile and adult survival") + ylab("")+
    scale_x_discrete(labels= c("gamma1"=expression(gamma[1]), "gamma2"=expression(gamma[2]), 
                               "phi1"=expression(phi[1]), "phi2"=expression(phi[2]), "s1a"=expression("s"["1a"]), "s2a"=expression("s"["2a"])))
  
  plot3 <- ggplot(borne_coexist[9:16,], aes(x=rownames(borne_coexist[9:16,]), y=borne_coexist[9:16,1])) + 
    geom_pointrange(ymin=(borne_coexist[9:16,1]+borne_coexist[9:16,2]), ymax=(borne_coexist[9:16,1]+borne_coexist[9:16,3]))+ 
    coord_cartesian(ylim = c(0, 0.4))+
    xlab("competition parameters") + ylab("")+
    scale_x_discrete(labels= c("alpha11"=expression(alpha[11]), "alpha12"=expression(alpha[12]), "alpha21"=expression(alpha[21]), "alpha22"=expression(alpha[22]),
                               "beta11"=expression(beta[11]), "beta12"=expression(beta[12]), "beta21"=expression(beta[21]), "beta22"=expression(beta[22])))
  
  a[[i]] <- plot_grid(plot1, plot2, plot3, ncol=3, nrow = 1, labels = LETTERS[i])
}

cairo_pdf("coexistence_domains_3sets.pdf", width = 10)
plot_grid(a[[1]], a[[2]], a[[3]], ncol=1, nrow = 3)
dev.off()

cairo_ps("coexistence_domains_3sets.eps", width = 10)
plot_grid(a[[1]], a[[2]], a[[3]], ncol=1, nrow = 3)
dev.off()
