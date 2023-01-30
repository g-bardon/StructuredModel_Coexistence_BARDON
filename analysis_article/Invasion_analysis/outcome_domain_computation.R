### GAEL BARDON ###
### 05/10/2020 ###
### Computation of the invasion criteria over alpha_ij and alpha_ji in the combined model with three situations suggested by beta coefficients 

library("gridExtra")
library(ggplot2)
library(cowplot)
library(Cairo)

## Function to compute the border of the outcome domains on alpha from invasion criteria
computation_border_combined <- function(pi1,pi2, gamma1,gamma2, phi1,phi2, sa1,sa2, alpha11, alpha22, beta11, beta12, beta21, beta22){
  n_a1 = (alpha11*(1-gamma1)*phi1 - alpha11 - beta11 + sqrt((-alpha11*(1-gamma1)*phi1 + alpha11 + beta11)^2-4*alpha11*beta11*(-(1-gamma1)*phi1-pi1*gamma1*phi1/(1-sa1)+1)))/(2*alpha11*beta11)
  alpha21_b = ((pi2*gamma2*phi2/(1-sa2))/((1+beta21*n_a1)-(1-gamma2)*phi2)-1)*1/n_a1
  n_a2 = (alpha22*(1-gamma2)*phi2 - alpha22 - beta22 + sqrt((-alpha22*(1-gamma2)*phi2 + alpha22 + beta22)^2-4*alpha22*beta22*(-(1-gamma2)*phi2-pi2*gamma2*phi2/(1-sa2)+1)))/(2*alpha22*beta22)
  alpha12_b = ((pi1*gamma1*phi1/(1-sa1))/((1+beta12*n_a2)-(1-gamma1)*phi1)-1)*1/n_a2
  return(c(alpha21_b,alpha12_b))
}

## Function to compute the border of the outcome domains from invasion criteria when beta=0 (From Fujiwara et al. but gives the same results if we use our method)
computation_border_only_alpha <- function(pi1,pi2, gamma1,gamma2, phi1,phi2, sa1,sa2, alpha11, alpha22){
  f1<-(1/(gamma1*phi1))*(1-sa1)*(1-phi1+gamma1*phi1)
  f2<-(1/(gamma2*phi2))*(1-sa2)*(1-phi2+gamma2*phi2)
  Rf1<-(pi1/f1)-1
  Rf2<-(pi2/f2)-1
  alpha12_b<-(Rf1/Rf2)*alpha22
  alpha21_b<-(Rf2/Rf1)*alpha11
  return(c(alpha21_b,alpha12_b))
} 

## Function to compute the border of the outcome domains from invasion criteria when beta=alpha 
computation_border_alpha_eq_beta <- function(pi1,pi2, gamma1,gamma2, phi1,phi2, sa1,sa2, alpha11, alpha22, beta11, beta22){
  n_a1 = (alpha11*(1-gamma1)*phi1 - alpha11 - beta11 + sqrt((-alpha11*(1-gamma1)*phi1 + alpha11 + beta11)^2-4*alpha11*beta11*(-(1-gamma1)*phi1-pi1*gamma1*phi1/(1-sa1)+1)))/(2*alpha11*beta11)
  n_a2 = (alpha22*(1-gamma2)*phi2 - alpha22 - beta22 + sqrt((-alpha22*(1-gamma2)*phi2 + alpha22 + beta22)^2-4*alpha22*beta22*(-(1-gamma2)*phi2-pi2*gamma2*phi2/(1-sa2)+1)))/(2*alpha22*beta22)
  C1= (1-gamma1)*phi1
  D1= pi1*gamma1*phi1/(1-sa1)
  C2= (1-gamma2)*phi2
  D2= pi2*gamma2*phi2/(1-sa2)
  alpha21_b = (sqrt(D2-C2)-1)*1/n_a1
  alpha12_b = (sqrt(D1-C1)-1)*1/n_a2
  return(c(alpha21_b,alpha12_b))
}


### Reference plot : inavsion criteria with only alpha (beta=0)
border_ref <- computation_border_only_alpha(30,25,0.7,0.8,0.5,0.4,0.5,0.6,0.1,0.1)
rect_ref <- data.frame(xmin = c(0,border_ref[1],0,border_ref[1]),
                      xmax = c(border_ref[1], 0.5,border_ref[1],0.5),
                      ymin = c(0 , 0, border_ref[2],border_ref[2]),
                      ymax = c(border_ref[2],border_ref[2],0.5,0.5),
                      grp = factor(c('Coexistence','Species 1 excludes species 2','Species 2 excludes species 1','Priority effect')))

a <- ggplot(data=rect_ref, aes(x=c(0,0.5), y=c(0,0.5)))+
  geom_rect(data = rect_ref,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, fill = grp)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('grey',"white","grey20", "grey50"))+
  xlab(expression(alpha[21]))+ylab(expression(alpha[12]))+ guides(fill=guide_legend(title="Outcome predicted by the invasion criteria"))+ggtitle("Competition on fertility only (\u03B1 = f(N))")+ 
  theme(plot.title = element_text(size = 12, face = "bold"))

### Reference beta =alpha

border_alpha_eq_beta <- computation_border_alpha_eq_beta(30,25,0.7,0.8,0.5,0.4,0.5,0.6,0.1,0.1,0.1,0.1)


rect_alpha_eq_beta <- data.frame(xmin = c(0,border_alpha_eq_beta[1],0,border_alpha_eq_beta[1]),
                      xmax = c(border_alpha_eq_beta[1], 0.5,border_alpha_eq_beta[1],0.5),
                      ymin = c(0 , 0, border_alpha_eq_beta[2],border_alpha_eq_beta[2]),
                      ymax = c(border_alpha_eq_beta[2],border_alpha_eq_beta[2],0.5,0.5),
                      grp = factor(c('Coexistence','Species 1 excludes species 2','Species 2 excludes species 1','Priority effect')))
b <- ggplot(data=rect_alpha_eq_beta, aes(x=c(0,0.5), y=c(0,0.5)))+
  geom_rect(data = rect_alpha_eq_beta,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, fill = grp)) +
  geom_segment(aes(x = 0, y = border_ref[2], xend = 0.5, yend = border_ref[2]), linetype="dashed", size=1)+
  geom_segment(aes(x = border_ref[1], y = 0, xend = border_ref[1], yend = 0.5), linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('grey',"white","grey20", "grey50"))+
  xlab(expression(alpha[21]))+ylab(expression(alpha[12]))+ guides(fill=guide_legend(title="Outcome predicted by the invasion criteria"))+ggtitle("Competition on both vital rates with \u03B2 = \u03B1")+ 
  theme(plot.title = element_text(size = 12, face = "bold"))+theme(legend.position = "none")


### beta suggests coexistence 
coexistence_border <- computation_border_combined(30,25,0.7,0.8,0.5,0.4,0.5,0.6,0.1,0.1,0.1, 0.06, 0.06, 0.1)

rect_coexistence <- data.frame(xmin = c(0,coexistence_border[1],0,coexistence_border[1]),
                                    xmax = c(coexistence_border[1], 0.5,coexistence_border[1],0.5),
                                    ymin = c(0 , 0, coexistence_border[2],coexistence_border[2]),
                                    ymax = c(coexistence_border[2],coexistence_border[2],0.5,0.5),
                                    grp = factor(c('Coexistence','Species 1 excludes species 2','Species 2 excludes species 1','Priority effect')))
c <- ggplot(data=rect_coexistence, aes(x=c(0,0.5), y=c(0,0.5)))+
  geom_rect(data = rect_coexistence,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, fill = grp)) +
  geom_segment(aes(x = 0, y = border_ref[2], xend = 0.5, yend = border_ref[2]), linetype="dashed", size=1)+
  geom_segment(aes(x = border_ref[1], y = 0, xend = border_ref[1], yend = 0.5), linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('grey',"white","grey20", "grey50"))+
  xlab(expression(alpha[21]))+ylab(expression(alpha[12]))+ guides(fill=guide_legend(title="Outcome predicted by the invasion criteria"))+ggtitle("\u03B2 indicates coexistence ")+ 
  theme(plot.title = element_text(size = 12, face = "bold"))+theme(legend.position = "none")+
  geom_point(aes(x=0.06,y=0.05), shape=8)

### Beta suggest exclusion sp 2

exclu_sp2_border <- computation_border_combined(30,25,0.7,0.8,0.5,0.4,0.5,0.6,0.1,0.1,0.1,0.125,0.01,0.1)

rect_exclu_sp2 <- data.frame(xmin = c(0,exclu_sp2_border[1],0,exclu_sp2_border[1]),
                               xmax = c(exclu_sp2_border[1], 0.5,exclu_sp2_border[1],0.5),
                               ymin = c(0 , 0, exclu_sp2_border[2],exclu_sp2_border[2]),
                               ymax = c(exclu_sp2_border[2],exclu_sp2_border[2],0.5,0.5),
                               grp = factor(c('Coexistence','Species 1 excludes species 2','Species 2 excludes species 1','Priority effect')))
d <- ggplot(data=rect_exclu_sp2, aes(x=c(0,0.5), y=c(0,0.5)))+
  geom_rect(data = rect_exclu_sp2,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, fill = grp)) +
  geom_segment(aes(x = 0, y = border_ref[2], xend = 0.5, yend = border_ref[2]), linetype="dashed", size=1)+
  geom_segment(aes(x = border_ref[1], y = 0, xend = border_ref[1], yend = 0.5), linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('grey',"white","grey20", "grey50"))+
  xlab(expression(alpha[21]))+ylab(expression(alpha[12]))+ guides(fill=guide_legend(title="Outcome predicted by the invasion criteria"))+ggtitle("\u03B2 indicates exclusion of species 1")+ 
  theme(plot.title = element_text(size = 12, face = "bold"))+theme(legend.position = "none")+
  geom_point(aes(x=0.112,y=0.02), shape=8)

#### BEta suggests priority effect
PE_border <- computation_border_combined(30,25,0.7,0.8,0.5,0.4,0.5,0.6,0.1,0.1,0.1,0.155,0.165,0.1)

rect_priorityeffect <- data.frame(xmin = c(0,PE_border[1],0,PE_border[1]),
                             xmax = c(PE_border[1], 0.5,PE_border[1],0.5),
                             ymin = c(0 , 0, PE_border[2],PE_border[2]),
                             ymax = c(PE_border[2],PE_border[2],0.5,0.5),
                             grp = factor(c('Coexistence','Species 1 excludes species 2','Species 2 excludes species 1','Priority effect')))
e <- ggplot(data=rect_priorityeffect, aes(x=c(0,0.5), y=c(0,0.5)))+
  geom_rect(data = rect_priorityeffect,aes(x = NULL,y = NULL,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, fill = grp)) +
  geom_segment(aes(x = 0, y = border_ref[2], xend = 0.5, yend = border_ref[2]), linetype="dashed", size=1)+
  geom_segment(aes(x = border_ref[1], y = 0, xend = border_ref[1], yend = 0.5), linetype="dashed", size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c('grey',"white","grey20", "grey50"))+
  xlab(expression(alpha[21]))+ylab(expression(alpha[12]))+ guides(fill=guide_legend(title="Outcome predicted by the invasion criteria"))+ggtitle("\u03B2 indicates a priority effect")+ 
  theme(plot.title = element_text(size = 12, face = "bold"))+theme(legend.position = "none")+
  geom_point(aes(x=0.035,y=0.043), shape=8)

#### Combined figurewith labels
f <- get_legend(a)
a <- a + theme(legend.position = "none")

cairo_pdf(file = "outcome_domains_full.pdf", width = 8)
plot_grid(a, b,c,d,e,f, labels=c("A", "B","C","D", "E"), ncol = 2, nrow = 3)
dev.off()


