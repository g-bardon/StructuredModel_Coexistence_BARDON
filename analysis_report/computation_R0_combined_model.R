### GAEL BARDON ###
### 21/07/2020 ###
### Computation of the invasion criteria in the combined model for both species in grid of parameter (which can be changed)



source("sensitivity_analysis_function.R")
library(ggplot2)
library('plot.matrix')


computation_R0 <- function(fert1,fert2, mat_rate1,mat_rate2, juv_surv1,juv_surv2, adult_surv1,adult_surv2, alpha11, alpha12, alpha21, alpha22, beta11, beta12, beta21, beta22){
  n_a1 = (alpha11*(1-mat_rate1)*juv_surv1 - alpha11 - beta11 + sqrt((-alpha11*(1-mat_rate1)*juv_surv1 + alpha11 + beta11)^2-4*alpha11*beta11*(-(1-mat_rate1)*juv_surv1-fert1*mat_rate1*juv_surv1/(1-adult_surv1)+1)))/(2*alpha11*beta11)
  R02 <- (1-mat_rate2)*juv_surv2/(1+beta21*n_a1)+(fert2*mat_rate2*juv_surv2/(1-adult_surv2))/((1+beta21*n_a1)*(1+alpha21*n_a1))
  n_a2 = (alpha22*(1-mat_rate2)*juv_surv2 - alpha22 - beta22 + sqrt((-alpha22*(1-mat_rate2)*juv_surv2 + alpha22 + beta22)^2-4*alpha22*beta22*(-(1-mat_rate2)*juv_surv2-fert2*mat_rate2*juv_surv2/(1-adult_surv2)+1)))/(2*alpha22*beta22)
  R01 <- (1-mat_rate1)*juv_surv1/(1+beta12*n_a2)+(fert1*mat_rate1*juv_surv1/(1-adult_surv1))/((1+beta12*n_a2)*(1+alpha12*n_a2))
  return(c(R01,R02))
}

l=seq(0.001,0.3,0.001)  ### Allow to define the range of the parameter
mat1=matrix(0,length(l),length(l), dimnames = list(l,l))
mat2=matrix(0,length(l),length(l), dimnames = list(l,l))
for (i in 1:length(l)){
  for (j in 1:length(l)){
    mat1[i,j]=computation_R0(5,5,0.5,0.5,0.5,0.5,0.5,0.5,0.1,l[i],0.05,0.1,0.1,l[j],0.05,0.1)[1]  ## only need to change the place of l[i] and l[j] to compute the invasion criteria over different parameter
    mat2[i,j]=computation_R0(5,5,0.5,0.5,0.5,0.5,0.5,0.5,0.1,l[i],0.05,0.1,0.1,l[j],0.05,0.1)[2]
  }
}

## Define the outcome predicted by both invasion criteria
coex_domain=matrix(0,length(l),length(l), dimnames = list(l,l))  
coex_domain[which(mat1 <=1 & mat2 <=1)] = "Priority Effect"
coex_domain[which(mat1 >1 & mat2 <=1)] = "Species 1 wins"
coex_domain[which(mat1 <=1 & mat2 >1)] = "Species 2 wins"
coex_domain[which(mat1 >1 & mat2 >1)] = "Coexistence"

coex_domain= as.data.frame(as.table(coex_domain))
coex_domain[,1] <- as.numeric(as.character(coex_domain[,1]))
coex_domain[,2] <- as.numeric(as.character(coex_domain[,2]))

### Plot of the domain
ggplot(coex_domain, aes(x=coex_domain[,1], y=coex_domain[,2]))+
  geom_raster(aes(fill=coex_domain[,3]), interpolate = FALSE)+
  xlim(0,0.3)+
  ylim(0,0.3)+
  xlab("alpha12")+ylab("alpha21")+ guides(fill=guide_legend(title="Outcome predicted by R0"))+ggtitle("Beta's suggest coexistence")

