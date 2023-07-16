rm(list=ls())

source("R/swarm_functions.R")

set.seed(123)
model_function <-function(params,ls){
  parammat <- matrix(params,nrow=3)
  l1_w <- parammat[1,]
  l1_b <- parammat[2,]
  lo_2 <- parammat[3,]

  modval<-t(pmax((l1_w%*%t(xdat))+l1_b,((l1_w%*%t(xdat))+l1_b)*0.03))%*%t(t(lo_2))
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

ls = 100 #Layer size
ndat=50
xdat=sort(runif(ndat,-10,10))
ydat = c(sin(0.5*xdat)/2)+0.6

lowlim=rep(-1,3*ls)
highlim=rep(1,3*ls)
#Test
swarm_state <- initialize_swarm(desired_values = ydat,param_len = 3*ls,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w = c())

for (itt in 1:50){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 1e-2)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w = c(),config = swarm.config(swarm.control = tempcontrol))
  matplot(xdat,swarm_state$pout,type="l",col="grey",ylim=c(0,1.2))
  points(xdat,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=1)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec),ls),col="red")
}

