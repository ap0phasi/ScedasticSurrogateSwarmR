rm(list=ls())

source("R/swarm_functions.R")

library(DiceKriging)


model_function <- function(params){
  modout <- rep(DiceKriging::hartman6(params),2)
  return(modout)
}

#Equality Constraints
eq <- function(params){
  constraints = c(
    # params[2]+params[5]-1.6
  )
  return(constraints)
}

#Inequality Constraints
# of form <= 0
ineq <- function(params){
  constraints = c(
    # params[1]-params[4]
  )
  return(constraints)
}

bestval = c(0.206,0.150011,0.476874,0.275332,0.311652,0.6573)
truemin = DiceKriging::hartman6(bestval)

ndat=1
ydat=c(-4,-4)
xplot = c(1,2)

lowlim=rep(0,6)
highlim=rep(1,6)
swarm_state <- initialize_swarm(desired_values = ydat,param_len = 6,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=21,num_cluster=3,deg=1))

for (itt in 1:20){
  tempcontrol = swarm.control(poly_w=0.5,stoch_w = 0.1)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=21,num_cluster=3,swarm.control = tempcontrol))
  print(min(apply(swarm_state$best_p,1,DiceKriging::hartman6)))
  matplot(t(swarm_state$best_p),type="l",col="grey")
  lines(bestval,col="red")
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=1)
}

#Compare to BOBYQA
result <- minqa::bobyqa(runif(6,0,1), DiceKriging::hartman6, lower = lowlim, upper = highlim)
lines(result$par,col="blue")
