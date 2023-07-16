rm(list=ls())

source("R/swarm_functions.R")

set.seed(128)

model_function <- function(params){
  parammat=matrix(params,nrow=3)
  centers=parammat[1,]
  weights=parammat[2,]
  widths=parammat[3,]*10
  modout=rep(0,length(xdat))
  for (ix in 1:length(centers)){
    modout=modout+weights[ix]*exp(-(widths[ix]*abs(xdat-centers[ix]))^2)
  }

  return(modout)
}

#Equality Constraints
eq <- function(params){
  constraints = c(
     params[2]+params[5]-1.6
  )
  return(constraints)
}

#Inequality Constraints
# of form <= 0
ineq <- function(params){
  constraints = c(
     params[1]-params[4]
  )
  return(constraints)
}

ndat=200
xdat=sort(runif(ndat,0,1))
params=c(0.15,0.8,0.5,0.8,0.8,0.5)
ydat=model_function(params)
ydat[xdat<0.3]=ydat[xdat<0.3]+rnorm(sum(xdat<0.3),0,0.1)
ydat = c(ydat,0,0)
xplot= c(xdat,seq(max(xdat),max(xdat)*1.1,length.out=2))

lowlim=rep(0,6)
highlim=rep(1.2,6)
#Test
swarm_state <- initialize_swarm(desired_values = ydat,param_len = 6,lowlim=lowlim,highlim=highlim,ineq_w = c(100),eq_w=c(100),config = swarm.config(num_cluster=3,deg=2))

for (itt in 1:100){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 100)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(100),eq_w=c(100),config = swarm.config(num_cluster=3,swarm.control = tempcontrol))
  matplot(xplot,swarm_state$pout,type="l",col="grey")
  points(xplot,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=2)
  lines(xplot,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec)),col="red")
}

bobya_function <-function(params){
  sum(abs(model_function(params)-ydat))
}
#Compare to BOBYA
result <- minqa::bobyqa(runif(6,0,1), bobya_function, lower = lowlim, upper = highlim)
lines(xdat,model_function(as.vector(result$par)),col="green")
