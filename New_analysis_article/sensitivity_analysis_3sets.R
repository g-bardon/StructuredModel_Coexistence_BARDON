library(ggplot2)
library(gridExtra)

### SENSITIVITY ANALYSIS REQUIRES PARAMETERS SET WITH THE FOLLOWING SHAPE : 
#(pi01, pi02, gamma1, gamma2, phi1, phi2, s1a, s2a, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22)
name_param = c("\u03c01", "\u03c02", "\u03b31", "\u03b32", "\u03c61", "\u03c62", "s1a", "s2a", "\u03b111", "\u03b112", "\u03b121", "\u03b122", "\u03b211", "\u03b212", "\u03b221", "\u03b222")

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
specie <- c(name_param, name_param, name_param)
data[,1] <- specie
data <- as.data.frame(data)
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))
# Grouped species 1
ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + 
  scale_x_discrete(limits = name_param) +
  ylab("Sensitivity") + xlab("") +
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Sensitivities of the parameters for the species 1")

# Grouped species 2
ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + 
  scale_x_discrete(limits = name_param) +
  ylab("Sensitivity") + xlab("") + 
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Sensitivities of the parameters for the species 2")


####

elas1 <- cbind(elas1, 1)
elas2 <- cbind(elas2, 2)
elas3 <- cbind(elas3, 3)
data <- rbind(elas1, elas2, elas3)
specie <- c(name_param, name_param, name_param)
data[,1] <- specie
data <- as.data.frame(data)
data[,2] <- as.numeric(as.character(data[,2]))
data[,3] <- as.numeric(as.character(data[,3]))
# Grouped species 1
ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + 
  scale_x_discrete(limits = name_param) +
  ylab("Elasticity") + xlab("") +
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Elasticities of the parameters for the species 1")

# Grouped species 2
ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
  geom_bar(position="dodge", stat="identity")  + 
  scale_x_discrete(limits = name_param) +
  ylab("Elasticity") + xlab("") + 
  guides(fill=guide_legend("Parameter set"))+
  labs(title = "Elasticities of the parameters for the species 2")

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
  rownames(borne_coexist) <- c("\u03c01", "\u03c02", "\u03b31", "\u03b32", "\u03c61", "\u03c62", "s1a", "s2a", "\u03b111", "\u03b112", "\u03b121", "\u03b122", "\u03b211", "\u03b212", "\u03b221", "\u03b222")
  plot1 <- ggplot(borne_coexist[1:2,], aes(x=rownames(borne_coexist[1:2,]), y=borne_coexist[1:2,1])) + 
    geom_pointrange(ymin=(borne_coexist[1:2,1]+borne_coexist[1:2,2]), ymax=(borne_coexist[1:2,1]+borne_coexist[1:2,3]))+ 
    coord_cartesian(ylim = c(0, 50))+
    xlab("Fertility") + ylab("Parameter value")
  
  plot2 <- ggplot(borne_coexist[3:8,], aes(x=rownames(borne_coexist[3:8,]), y=borne_coexist[3:8,1])) + 
    geom_pointrange(ymin=(borne_coexist[3:8,1]+borne_coexist[3:8,2]), ymax=(borne_coexist[3:8,1]+borne_coexist[3:8,3]))+ 
    coord_cartesian(ylim = c(0, 1))+
    xlab("Maturation rate, juvenile and adult survival") + ylab("")
  
  plot3 <- ggplot(borne_coexist[9:16,], aes(x=rownames(borne_coexist[9:16,]), y=borne_coexist[9:16,1])) + 
    geom_pointrange(ymin=(borne_coexist[9:16,1]+borne_coexist[9:16,2]), ymax=(borne_coexist[9:16,1]+borne_coexist[9:16,3]))+ 
    coord_cartesian(ylim = c(0, 0.4))+
    xlab("competition parameters") + ylab("")
  
  a[[i]] <- plot_grid(plot1, plot2, plot3, ncol=3, nrow = 1, labels = LETTERS[i])
}

plot_grid(a[[1]], a[[2]], a[[3]], ncol=1, nrow = 3)
