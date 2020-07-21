### GAEL BARDON ###
### 21/07/2020 ###
### Representation of the intersection of conic sections that come from the equation allowing to compute the densities at equilbirum for the structured combined model of competition.

alpha11 <- runif(1, 0.5,2)
beta12 <- runif(1, 0.2,1)
alpha12 <- runif(1, 0.2,1)
beta11 <- runif(1, 0.5,2)
c1 <- runif(1,0.2,0.9)
d1 <- runif(1,1,5)

alpha22 <- runif(1, 0.5,2)
beta21 <- runif(1, 0.2,1)
alpha21 <- runif(1, 0.2,1)
beta22 <- runif(1, 0.5,2)
c2 <- runif(1,0.2,0.9)
d2 <- runif(1,1,5)


## Function f1 and f2 are the conics sections and the intersection points are the solutions of the system
f1 <- function(x,y) alpha11*beta11*x^2+alpha12*beta12*y^2+(alpha11*beta12+alpha12*beta11)*x*y+(beta11+alpha11-c1*alpha11)*x+(beta12+alpha12-c1*alpha12)*y+(1-c1-d1)
x <- y <- seq(-5,5,length=100)
z <- outer(x,y,f1)
contour(
  x=x, y=x, z=z, 
  levels=0, las=1, drawlabels=FALSE, lwd=3,panel.first = grid(), col="green", xlab="x", ylab='y'
)
lines(c(0,0),c(5,-5))
lines(c(-5,5),c(0,0))

f2 <- function(x,y) alpha22*beta22*y^2+alpha21*beta21*x^2+(alpha22*beta21+alpha21*beta22)*x*y+(beta22+alpha22-c2*alpha22)*y+(beta21+alpha21-c2*alpha21)*x+(1-c2-d2)
x <- y <- seq(-5,5,length=100)
z <- outer(x,y,f2)
contour(
  x=x, y=x, z=z, 
  levels=0, las=1, drawlabels=FALSE, lwd=3, add=TRUE,panel.first = grid(), col="red"
)
legend("topright", legend=c("Conic section for species 1", "Conic section for species 2"), col=c("green", "red"), lty=1)

