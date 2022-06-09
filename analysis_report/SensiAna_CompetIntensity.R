library(ggplot2)
library(gridExtra)

### SENSITIVITY ANALYSIS REQUIRES PARAMETERS SET WITH THE FOLLOWING SHAPE : 
#(pi01, pi02, gamma1, gamma2, phi1, phi2, s1a, s2a, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22)

name_param = c("\u03c01", "\u03c02", "\u03b31", "\u03b32", "\u03c61", "\u03c62", "s1a", "s2a", "\u03b111", "\u03b112", "\u03b121", "\u03b122", "\u03b211", "\u03b212", "\u03b221", "\u03b222")

tab=list()
tab[[1]] <- rbind(c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.5,0.1,0.1,0.5,0.5,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.1,0.1,0.02,0.5,0.5,0.02,0.1,0.1),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.1,0.1,0.02,0.5,0.1,0.1,0.02,0.5))

tab[[2]] <- rbind(c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.2,0.1,0.1,0.2,0.2,0.1,0.1,0.2),
                   c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.04,0.1,0.02,0.2,0.2,0.02,0.1,0.04),
                   c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.04,0.1,0.02,0.2, 0.04,0.1,0.02,0.2))

tab[[3]] <- rbind(c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.11,0.1,0.1,0.11,0.11,0.1,0.1,0.11),
                   c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.022,0.1,0.02,0.11,0.11,0.02,0.1,0.022),
                   c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.022,0.1,0.02,0.11,0.022,0.1,0.02,0.11))

source("sensitivity_analysis_function.R")
ratio = c(5,2,1.1)

for (i in 1:3){
  tab_param = tab[[i]]
  elas1 <- sensitivity(tab_param[1,])$Elasticity
  elas2 <- sensitivity(tab_param[2,])$Elasticity
  elas3 <- sensitivity(tab_param[3,])$Elasticity
  
  sensi1 <- sensitivity(tab_param[1,])$Sensitivity
  sensi2 <- sensitivity(tab_param[2,])$Sensitivity
  sensi3 <- sensitivity(tab_param[3,])$Sensitivity
  

  sensi1 <- cbind(sensi1, 3*(i-1)+1)
  sensi2 <- cbind(sensi2, 3*(i-1)+2)
  sensi3 <- cbind(sensi3, 3*(i-1)+3)
  data <- rbind(sensi1, sensi2, sensi3)
  specie <- c(name_param, name_param, name_param)
  data[,1] <- specie
  data <- as.data.frame(data)
  data[,2] <- as.numeric(as.character(data[,2]))
  data[,3] <- as.numeric(as.character(data[,3]))
  # Grouped species 1
  print(ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
    geom_bar(position="dodge", stat="identity")  + 
    scale_x_discrete(limits = name_param) +
    ylab("Sensitivity") + xlab("") +
    guides(fill=guide_legend("Parameter set"))+
    labs(title = paste("Sensitivities of the parameters for the species 1 with a ratio of intra/inter competition coefficient equal to", ratio[i], sep="")))
  
  # Grouped species 2
  print(ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
    geom_bar(position="dodge", stat="identity")  + 
    scale_x_discrete(limits = name_param) +
    ylab("Sensitivity") + xlab("") + 
    guides(fill=guide_legend("Parameter set"))+
    labs(title = paste("Sensitivities of the parameters for the species 2 with a ratio of intra/inter competition coefficient equal to", ratio[i], sep="")))
  
  
  ####
  
  elas1 <- cbind(elas1, 3*(i-1)+1)
  elas2 <- cbind(elas2, 3*(i-1)+2)
  elas3 <- cbind(elas3, 3*(i-1)+3)
  data <- rbind(elas1, elas2, elas3)
  specie <- c(name_param, name_param, name_param)
  data[,1] <- specie
  data <- as.data.frame(data)
  data[,2] <- as.numeric(as.character(data[,2]))
  data[,3] <- as.numeric(as.character(data[,3]))
  # Grouped species 1
  print(ggplot(data, aes(fill=data[,4], y=data[,2], x=data[,1])) + 
    geom_bar(position="dodge", stat="identity")  + 
    scale_x_discrete(limits = name_param) +
    ylab("Elasticity") + xlab("") +
    guides(fill=guide_legend("Parameter set"))+
    labs(title = paste("Elasticities of the parameters for the species 1 with a ratio of intra/inter competition coefficient equal to", ratio[i])))
  
  # Grouped species 2
  print(ggplot(data, aes(fill=data[,4], y=data[,3], x=data[,1])) + 
    geom_bar(position="dodge", stat="identity")  + 
    scale_x_discrete(limits = name_param) +
    ylab("Elasticity") + xlab("") + 
    guides(fill=guide_legend("Parameter set"))+
    labs(title = paste("Elasticities of the parameters for the species 2 with a ratio of intra/inter competition coefficient equal to", ratio[i])))
  
  
  print(sensitivity(tab_param[1,])$Determinant)
  print(sensitivity(tab_param[2,])$Determinant)
  print(sensitivity(tab_param[3,])$Determinant)
  
}


