# delete everything
rm(list=ls(all=TRUE))

# import library for solving differential equations
require(deSolve)

shortEndTime <- 120

thisParams <- list(
  rho = 6000,
  # Scenario B (balanced primary and secondary)
  beta = 10^-3.92,
  epsilon = 10^-0.3,
  nu = 1,
  eta = 1/3,
  gamma = 1/9,
  alpha = 1/2,
  omega = 2/9,
  sigma = 1/2,
  mu = 1,
  # Conversion factor for external inoc is the number of infected leaves at eqm in the single-grower model
  incConvFact = 5667.40326042357,
  delta = 0.2,
  propControl = 0.5 # called pi in the paper
)

#
# right hand side of area wide model
#
# Note there are now two populations; _C (controllers) _NC (non-controllers)
#
# The fraction of each is propControl
#
modelFuncII <- function(t, y, p) 
{
  # state variables for the controllers
  S_F_C <- y[1]
  A_F_C <- y[2]
  S_M_C <- y[3]
  A_M_C <- y[4]
  V_M_C <- y[5]
  E_C   <- y[6]
  I_C   <- y[7]

  # state variables for the non-controllers
  S_F_NC <- y[8]
  A_F_NC <- y[9]
  S_M_NC <- y[10]
  A_M_NC <- y[11]
  V_M_NC <- y[12]
  E_NC   <- y[13]
  I_NC   <- y[14]
  
  dydt <- c()
  
  priRate_C <- (p$epsilon/p$incConvFact) * (p$propControl * I_C + (1-p$propControl) * I_NC)
  secRate_C <- p$beta * I_C

  priRate_NC <- (p$epsilon/p$incConvFact) * (p$propControl * I_C + (1-p$propControl) * I_NC)
  secRate_NC <- p$beta * I_NC

  dydt[1] <- p$rho - priRate_C * p$nu * S_F_C - secRate_C * p$nu * S_F_C - p$eta * S_F_C
  dydt[2] <- priRate_C * p$nu * S_F_C + secRate_C * p$nu * S_F_C - p$eta * A_F_C
  dydt[3] <- p$eta * S_F_C - priRate_C * S_M_C - secRate_C * S_M_C - p$gamma * S_M_C
  dydt[4] <- priRate_C * S_M_C + secRate_C * S_M_C - p$gamma * A_M_C - p$alpha * A_M_C
  dydt[5] <- p$alpha * A_M_C + p$eta * A_F_C - p$gamma * V_M_C - p$omega * V_M_C
  dydt[6] <- p$gamma * V_M_C + p$omega * V_M_C - p$sigma * E_C - p$delta * E_C
  dydt[7] <- p$sigma * E_C - p$mu * I_C - p$delta * I_C

  dydt[8] <- p$rho - priRate_NC * p$nu * S_F_NC - secRate_NC * p$nu * S_F_NC - p$eta * S_F_NC
  dydt[9] <- priRate_NC * p$nu * S_F_NC + secRate_NC * p$nu * S_F_NC - p$eta * A_F_NC
  dydt[10] <- p$eta * S_F_NC - priRate_NC * S_M_NC - secRate_NC * S_M_NC - p$gamma * S_M_NC
  dydt[11] <- priRate_NC * S_M_NC + secRate_NC * S_M_NC - p$gamma * A_M_NC - p$alpha * A_M_NC
  dydt[12] <- p$alpha * A_M_NC + p$eta * A_F_NC - p$gamma * V_M_NC - p$omega * V_M_NC
  dydt[13] <- p$gamma * V_M_NC + p$omega * V_M_NC - p$sigma * E_NC  # (note non-controllers don't remove litter)  - p$delta * E_NC
  dydt[14] <- p$sigma * E_NC - p$mu * I_NC # (note non-controllers don't remove litter) - p$delta * I_NC
    
  return(list(dydt))  
}

initCond <- rep(c(thisParams$rho/thisParams$eta,0,thisParams$rho/thisParams$gamma,0,0,1,0),2)

output <- lsoda(initCond, seq(from = 0, to = shortEndTime, length.out = 1000), modelFuncII, thisParams)

plot(output[, 1], cex.axis = 1.25, cex.lab  = 1.25, output[, 2], lwd=3, col = "black", type = "l", ylim = c(0,30000), xlab="Time (months)", ylab="Number of leaves", main="Controllers (Scenario B)") # S
lines(output[, 1], output[, 3],lwd=3, pch = 2,  lty = 2, col = "black") # I
lines(output[, 1], output[, 4], lwd=3, lty = 3, col = "black") # X
lines(output[, 1], output[, 5], lwd=3, lty = 4, col = "black") # Z
lines(output[, 1], output[, 6], lwd=3, lty = 5, col = "black") # Z

legend("topright", lty = c(1, 2, 3, 4, 5), lwd=3,
       col = c("black"), cex=1.25,
       legend = c(expression("Susceptible flush (S"[F]^{"C"}*")"), expression("Asymptomatic flush (A"[F]^{"C"}*")"),  expression("Susceptible mature (S"[M]^{"C"}*")"), expression("Asymptomatic mature (A"[M]^{"C"}*")"), expression("Symptomatic mature (V"[M]^{"C"}*")")))

ncOffset <- 7
plot(output[, 1], cex.axis = 1.25, cex.lab  = 1.25, output[, 2+ncOffset], lwd=3, col = "black", type = "l", ylim = c(0,30000), xlab="Time (months)", ylab="Number of leaves", main="Non-controllers (Scenario B)") # S
lines(output[, 1], output[, 3+ncOffset],lwd=3, pch = 2,  lty = 2, col = "black") # I
lines(output[, 1], output[, 4+ncOffset], lwd=3, lty = 3, col = "black") # X
lines(output[, 1], output[, 5+ncOffset], lwd=3, lty = 4, col = "black") # Z
lines(output[, 1], output[, 6+ncOffset], lwd=3, lty = 5, col = "black") # Z

legend("topright", lty = c(1, 2, 3, 4, 5), lwd=3,
       col = c("black"), cex=1.25,
       legend = c(expression("Susceptible flush (S"[F]^{"NC"}*")"), expression("Asymptomatic flush (A"[F]^{"NC"}*")"),  expression("Susceptible mature (S"[M]^{"NC"}*")"), expression("Asymptomatic mature (A"[M]^{"NC"}*")"), expression("Symptomatic mature (V"[M]^{"NC"}*")")))
