rm(list=ls())

source("R/swarm_functions.R")
source("R/search_functions.R")

#set.seed(128)

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

search_state <- initialize_search(param_len = 6,lowlim=lowlim,highlim=highlim,config = search.config(deg=1))

genstate = F
for (itt in 1:10){
  tempconfig = search.config(gen = genstate,deg=1,search_samples = 6*2,search_mag=0.3)
  search_state <- step_search(desired_values = ydat,search_state,ineq_w = c(),eq_w=c(),config = tempconfig)

  modout = model_function(as.vector(search_state$poly_recs[1,]))
  surrogate_out <- surrogate_model(as.vector(search_state$poly_recs[1,]),search_state$polyouts,search_state$centersaves,deg=1)
  if (mean(abs((modout-surrogate_out)/modout))>0.5){
    genstate = T
  }else{
    genstate = F
  }
  matplot(t(search_state$pos.x),type="l",col="grey")
  points(bestval)
  print(min(apply(search_state$pos.x,1,DiceKriging::hartman6)))
}
#usevec = (dim(search_state$xpins)[1]-200):dim(search_state$xpins)[1]
usevec = sample(1:dim(search_state$xpins)[1],12)

print(paste("Number of function evals: ",dim(search_state$xpins)[1]))

#Swap into particle swarm
swarm_state <- convert_search_to_swarm(search_state = search_state,usevec = usevec,desired_values = ydat,param_len = 6,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=length(usevec),num_cluster=1,deg=1))
#swarm_state <- initialize_swarm(desired_values = ydat,param_len = 3*ls,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w = c(),config = swarm.config(swarm_size = ls*2))

for (itt in 1:30){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 0.01)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=dim(swarm_state$x.p)[1],num_cluster=1,swarm.control = tempcontrol))
  print(min(apply(swarm_state$best_p,1,DiceKriging::hartman6)))
  matplot(t(swarm_state$best_p),type="l",col="grey")
}

print(paste("Number of function evals: ",dim(swarm_state$allxp)[1]))

print("-----BOBYQA-----")
bobya_function <-function(params){
  sum(abs(model_function(params)-ydat))
}
#Compare to BOBYA
result <- minqa::bobyqa(runif(6,0,1), bobya_function, lower = lowlim, upper = highlim)
lines(result$par,col="blue")
print(result)
points(bestval)
