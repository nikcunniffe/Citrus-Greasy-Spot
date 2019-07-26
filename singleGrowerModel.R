#
# Delete all variables from workspace
#
rm(list=ls(all=TRUE))

#
# Library for solving differential equations
#
require(deSolve)

#
# Set working directory
#
setwd(<YOUR DIRECTORY HERE>)

#
# This script replicates Fig. 10B and writes the results to a file
#

endTime <- 36

thisParams <- list(
  rho = 6000,
  beta = 10^-3.92,
  epsilon = 10^-0.3,
  nu = 1,
  eta = 1/3,
  gamma = 1/9,
  alpha = 1/2,
  omega = 2/9,
  sigma = 1/2,
  mu = 1,
  delta = 0
)

# right hand side of full differential equation
modelFunc <- function(t, y, p) 
{
  S_F <- y[1]
  A_F <- y[2]
  S_M <- y[3]
  A_M <- y[4]
  V_M <- y[5]
  E   <- y[6]
  I   <- y[7]

  dydt <- c()
  
  dydt[1] <- p$rho - p$epsilon * p$nu * S_F - p$beta * p$nu * I * S_F - p$eta * S_F
  dydt[2] <- p$epsilon * p$nu * S_F + p$beta * p$nu * I * S_F - p$eta * A_F
  dydt[3] <- p$eta * S_F - p$epsilon * S_M - p$beta * S_M * I - p$gamma * S_M
  dydt[4] <- p$epsilon * S_M + p$beta * S_M * I - p$gamma * A_M - p$alpha * A_M
  dydt[5] <- p$alpha * A_M + p$eta * A_F - p$gamma * V_M - p$omega * V_M
  dydt[6] <- p$gamma * V_M + p$omega * V_M - p$sigma * E - p$delta * E
  dydt[7] <- p$sigma * E - p$mu * I - p$delta * I
  
  return(list(dydt))  
}

png(file = "fig10B.png",width=8,height=8, units='in',res=300)

par(mar=c(5,7,2,2), las = 1, xpd = FALSE)
output <- lsoda(c(thisParams$rho/thisParams$eta,0,thisParams$rho/thisParams$gamma,0,0,1,0), seq(from = 0, to = 36, length.out = 1000), modelFunc, thisParams)

plot(output[, 1], cex.axis = 1.5, cex.lab  = 1.5, output[, 2], lwd=4, col = "#B2DF8A", type = "l",ylim = c(0,30000), xlab="Time (months)", ylab="", main="Baseline behaviour (Scenario B: Balanced infection)", cex.main = 1.5) # S

title(ylab="Number of leaves",line=5.5,cex.lab=1.5)

lines(output[, 1], output[, 3],lwd=4, pch = 2,  lty = 2, col = "#A6CEE3") # I
lines(output[, 1], output[, 4], lwd=4, lty = 3, col = "#33A02C") # X
lines(output[, 1], output[, 5], lwd=4, lty = 4, col = "#1F78B4") # Z
lines(output[, 1], output[, 6], lwd=4, lty = 5, col = "#984EA3") # Z

legend("topright", lty = c(1, 2, 3, 4, 5), lwd=3,
       col = c("#B2DF8A","#A6CEE3","#33A02C","#1F78B4","#984EA3"), cex=1.2,
       legend = c(expression("Susceptible flush (S"[F]*")"), expression("Asymptomatic flush (A"[F]*")"),  expression("Susceptible mature (S"[M]*")"), expression("Asymptomatic mature (A"[M]*")"), expression("Symptomatic mature (V"[M]*")")))

dev.off()

