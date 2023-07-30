rm(list=ls())

source("R/swarm_functions.R")
source("R/search_functions.R")

model_function <-function(params,ls){
  modval<-xdat*params[1]+params[2]
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
swarm_state <- initialize_search(param_len = 2,lowlim=lowlim,highlim=highlim,config = search.config())

genstate = F
for (itt in 1:10){
  tempconfig = search.config(gen = genstate)
  swarm_state <- step_search(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = tempconfig)
  plot(swarm_state$xpins)
  print(swarm_state$pos.x)

  modout = model_function(as.vector(swarm_state$poly_recs[1,]),ls)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_recs[1,]),swarm_state$polyouts,swarm_state$centersaves,deg=2)
  if (mean(abs((modout-surrogate_out)/modout))>0.1){
    genstate = T
  }else{
    genstate = F
  }
  # plot(xdat,ydat)
  # lines(xdat,surrogate_out,col="blue")
  # lines(xdat,model_function(as.vector(swarm_state$poly_recs[1,]),ls),col="red")
}

