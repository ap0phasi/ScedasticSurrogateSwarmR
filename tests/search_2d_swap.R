rm(list=ls())

source("R/swarm_functions.R")
source("R/search_functions.R")

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
search_state <- initialize_search(param_len = 2,lowlim=lowlim,highlim=highlim,config = search.config())

genstate = F
for (itt in 1:10){
  tempconfig = search.config(gen = genstate,revert_best = T)
  search_state <- step_search(desired_values = ydat,search_state,ineq_w = c(),eq_w=c(),config = tempconfig)
  plot(search_state$xpins)

  modout = model_function(as.vector(search_state$poly_recs[1,]),ls)
  surrogate_out <- surrogate_model(as.vector(search_state$poly_recs[1,]),search_state$polyouts,search_state$centersaves,deg=2)
  if (mean(abs((modout-surrogate_out)/modout))>0.5){
    genstate = T
  }else{
    genstate = F
  }
  matplot(xdat,search_state$current_pout,type="l",col="grey")
  points(xdat,ydat)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(search_state$poly_recs[1,]),ls),col="red")
}
print(paste("Number of function evals: ",dim(search_state$xpins)[1]))

usevec = (dim(search_state$xpins)[1]-50):dim(search_state$xpins)[1]

#Swap into particle swarm
swarm_state <- convert_search_to_swarm(search_state = search_state,usevec = usevec,desired_values = ydat,param_len = 2,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=length(usevec),num_cluster=3,deg=2))

for (itt in 1:10){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 0.01)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=length(usevec),num_cluster=3,swarm.control = tempcontrol))
  plot(swarm_state$x.p)
  matplot(xdat,swarm_state$pout,type="l",col="grey")
  points(xdat,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=2)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec),ls),col="red")
}
print(paste("Number of function evals: ",dim(swarm_state$allxp)[1]))

print("-----BOBYQA-----")
bobya_function <-function(params){
  sum(abs(model_function(params)-ydat))
}
#Compare to BOBYA
result <- minqa::bobyqa(runif(2,0,1), bobya_function, lower = lowlim, upper = highlim)
lines(xdat,model_function(as.vector(result$par)),col="green")
print(result)
