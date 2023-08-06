rm(list=ls())

source("R/swarm_functions.R")
source("R/search_functions.R")

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

ndat=200
xdat=sort(runif(ndat,0,1))
params=c(0.15,0.8,0.5,0.8,0.8,0.5)
ydat=model_function(params)
ydat[xdat<0.3]=ydat[xdat<0.3]+rnorm(sum(xdat<0.3),0,0.1)
xplot= c(xdat,seq(max(xdat),max(xdat)*1.1,length.out=0))

lowlim=rep(0,6)
highlim=rep(1.2,6)
#Test
search_state <- initialize_search(param_len = 6,lowlim=lowlim,highlim=highlim,config = search.config(deg=1))

genstate = F
for (itt in 1:7){
  tempconfig = search.config(gen = genstate,deg=1,search_samples = 6*2,search_mag=0.5)
  search_state <- step_search(desired_values = ydat,search_state,ineq_w = c(),eq_w=c(),config = tempconfig)

  modout = model_function(as.vector(search_state$poly_recs[1,]))
  surrogate_out <- surrogate_model(as.vector(search_state$poly_recs[1,]),search_state$polyouts,search_state$centersaves,deg=1)
  if (mean(abs((modout-surrogate_out)/modout))>0.5){
    genstate = T
  }else{
    genstate = F
  }
  matplot(xdat,search_state$current_pout,type="l",col="grey",ylim = c(0,1.2))
  points(xdat,ydat)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(search_state$poly_recs[1,])),col="red")
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
  matplot(xdat,swarm_state$pout,type="l",col="grey",ylim = c(0,1.2))
  points(xdat,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=1)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec)),col="red")
}

print(paste("Number of function evals: ",dim(swarm_state$allxp)[1]))

print("-----BOBYQA-----")
bobya_function <-function(params){
  sum(abs(model_function(params)-ydat))
}
#Compare to BOBYA
result <- minqa::bobyqa(runif(6,0,1.2), bobya_function, lower = lowlim, upper = highlim)
lines(xdat,model_function(as.vector(result$par)),col="green")
print(result)
