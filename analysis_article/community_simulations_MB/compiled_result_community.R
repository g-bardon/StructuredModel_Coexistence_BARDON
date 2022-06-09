library("ggplot2")
source("simul_functions.R")

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

var_list<-c("low", "average", "high")
S_list <- c("5", "10", "40")
c3 <- LETTERS[1:9]
gen=3000
plot_list <- list()
z=1
par(mfrow=c(3,3))
for (i in var_list){
  for (j in S_list){
    S=as.numeric(j)
    ini <- read.csv(paste(j,"_",i,"/initialization.csv", sep = ""))
    param <- read.csv(paste(j,"_",i,"/param_vital_rates.csv", sep=""))
    alpha <- as.matrix(read.table(paste(j,"_",i,"/matrix_alpha.csv", sep="")))
    beta <- as.matrix(read.table(paste(j,"_",i,"/matrix_beta.csv", sep="")))
    
    x_ref <- algo_simul(gen, S, ini, param[,1],param[,2],param[,3],param[,4],alpha,beta)
    color= rep("white", S)
    color[x_ref] = "gray"
    x <- simul_permut(gen, S, ini, param[,1] ,param[,2],param[,3],param[,4], alpha, beta)
    x<- c(x,x_ref)
    hist(x, xlim=c(0,S+0.5),breaks=0:S+0.5, col=color, xlab = "Number of species at equilibrium", main = paste(S, "species,",i,"variance"))
    fig_label(c3[z], cex=2)
    z=z+1
  }
}
full_plot = recordPlot()
#plot.new()

cairo_pdf("9hist_opp_hier_only_MB.pdf", height = 6)
full_plot
dev.off()

cairo_ps("9hist_opp_hier_only_MB.eps", height = 6)
full_plot
dev.off()
