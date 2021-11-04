library(bezier)

my.points <- matrix(c(0,1,2,2,4,1.5,5,3), byrow = TRUE, ncol=2 )




bezier.interpolate <- function(m, cp.dist.type="equal", cp.dist=0.5, npoints=100) {
  
  if(cp.dist.type == "equal") {
    cp1 <- m
    cp1[,1] <- cp1[,1] + cp.dist
    cp1 <- cp1[-nrow(cp1),]
    
    cp2 <- m
    cp2[,1] <- cp2[,1] - cp.dist
    cp2 <- cp2[-1,]
  } else {
    x.dist <- c(m[,1],0)-c(0, m[,1])
    x.dist <- x.dist[2:(nrow(m))]
    
    cp1 <- m
    cp1[,1] <- cp1[,1] + c(x.dist*cp.dist,0)
    cp1 <- cp1[-nrow(cp1),]
    
    cp2 <- m
    cp2[,1] <- cp2[,1] - c(0,x.dist*cp.dist)
    cp2 <- cp2[-1,]
  }
  
  l <- list(my.points, cp1, cp2)
  bz.matrix <- do.call(rbind, l)[order(sequence(sapply(l, nrow))), ]
  
  t <- seq(0, nrow(m)-1, length=100)
  
  return(bezier(t=t, p=bz.matrix, deg = 3))
}



plot(my.points, type = "p", xlim = c(0,5), ylim=c(0,4))
bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
lines(bz.points, type = "l")

my.points = m


plot(my.points, type = "p", xlim = c(0,5), ylim=c(0,4))
bz.points <- bezier.interpolate(my.points, cp.dist.type = "equal", cp.dist=0.5)
lines(bz.points, type = "l")


