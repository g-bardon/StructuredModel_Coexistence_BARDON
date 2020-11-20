library(ggplot2)
library(gridExtra)

### SENSITIVITY ANALYSIS REQUIRES PARAMETERS SET WITH THE FOLLOWING SHAPE : 
#(pi01, pi02, gamma1, gamma2, phi1, phi2, s1a, s2a, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22)

name_param = c("\u03c01", "\u03c02", "\u03b31", "\u03b32", "\u03c61", "\u03c62", "s1a", "s2a", "\u03b111", "\u03b112", "\u03b121", "\u03b122", "\u03b211", "\u03b212", "\u03b221", "\u03b222")

tab_param <- rbind(c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 0.5,0.1,0.1,0.5,0.5,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 1,0.1,0.1,0.5,1,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 2,0.1,0.1,0.5,2,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 5,0.1,0.1,0.5,5,0.1,0.1,0.5),
                  c(30,30,0.8,0.8,0.5,0.5,0.5,0.5, 10,0.1,0.1,0.5,10,0.1,0.1,0.5))


source("sensitivity_analysis_function.R")
ratio = c(5,2,1.1)
sensi = list()
elas = list()
data_sensi = c()
data_elas = c()

for (i in 1:5){
  print(colSums(equi_simul(tab_param[i,])))
  
  elas[[i]] <- sensitivity(tab_param[i,])$Elasticity
  sensi[[i]] <- sensitivity(tab_param[i,])$Sensitivity
  
  sensi[[i]][,1] <- name_param 
  sensi[[i]] <- cbind(sensi[[i]], i)
  data_sensi <- rbind(data_sensi, sensi[[i]])
  
  data_sensi <- as.data.frame(data_sensi)
  data_sensi[,2] <- as.numeric(as.character(data_sensi[,2]))
  data_sensi[,3] <- as.numeric(as.character(data_sensi[,3]))
  
  ###
  
  elas[[i]][,1] <- name_param 
  elas[[i]] <- cbind(elas[[i]], i)
  data_elas <- rbind(data_elas, elas[[i]])
  
  data_elas <- as.data.frame(data_elas)
  data_elas[,2] <- as.numeric(as.character(data_elas[,2]))
  data_elas[,3] <- as.numeric(as.character(data_elas[,3]))
  
  print(sensitivity(tab_param[i,])$Determinant)
  print(sensitivity(tab_param[i,])$Mij_matrix)
  print(solve(sensitivity(tab_param[i,])$Mij_matrix))
}



  # Grouped species 1
  p1 <- ggplot(data_sensi, aes(fill=data_sensi[,4], y=data_sensi[,2], x=data_sensi[,1])) +
          geom_bar(position="dodge", stat="identity")  +
          scale_x_discrete(limits = name_param) +
          ylab("Sensitivity") + xlab("") +
          guides(fill=FALSE)+
          labs(title = paste("Sensitivity for the species 1"))

  # Grouped species 2
  p2 <- ggplot(data_sensi, aes(fill=data_sensi[,4], y=data_sensi[,3], x=data_sensi[,1])) +
          geom_bar(position="dodge", stat="identity")  +
          scale_x_discrete(limits = name_param) +
          ylab("Sensitivity") + xlab("") +
          guides(fill=guide_legend("Parameter set"))+
          labs(title = paste("Sensitivity for the species 2"))


  ####

  
  # Grouped species 1
  p3 <- ggplot(data_elas, aes(fill=data_elas[,4], y=data_elas[,2], x=data_elas[,1])) +
          geom_bar(position="dodge", stat="identity")  +
          scale_x_discrete(limits = name_param) +
          ylab("Elasticity") + xlab("") +
          guides(fill=FALSE)+
          labs(title = paste("Elasticity for the species 1"))
  
  # Grouped species 2
  p4 <- ggplot(data_elas, aes(fill=data_elas[,4], y=data_elas[,3], x=data_elas[,1])) +
          geom_bar(position="dodge", stat="identity")  +
          scale_x_discrete(limits = name_param) +
          ylab("Elasticity") + xlab("") +
          guides(fill=guide_legend("Parameter set"))+
          labs(title = paste("Elasticity for the species 2"))
  # 
  
grid.arrange(p1,p2,p3,p4, widths=c(45,55))


