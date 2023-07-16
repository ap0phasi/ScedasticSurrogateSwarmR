rm(list=ls())

source("R/swarm_functions.R")

set.seed(128)
model_function <-function(params,ls){
  modval<-params[1]*cos(xdat*params[2])+params[2]*sin(xdat*params[1])
  return(modval)
}

#Equality Constraints
eq <- function(params){
  constraints = c(
  )
  return(constraints)
}

#Inequality Constraints
# of form <= 0
ineq <- function(params){
  constraints = c(
  )
  return(constraints)
}

xdat=seq(1,100,0.5)
params=c(0.32,0.4)
ydat=modf(params)

lowlim=rep(0,2)
highlim=rep(1,2)
#Test
swarm_state <- initialize_swarm(desired_values = ydat,param_len = 2,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w=c(),config = swarm.config(num_cluster=3,deg=2))

for (itt in 1:10){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 0.01)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = swarm.config(num_cluster=3,swarm.control = tempcontrol))
  matplot(xdat,swarm_state$pout,type="l",col="grey")
  points(xdat,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=2)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec),ls),col="red")
}

