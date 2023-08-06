rm(list=ls())

source("R/swarm_functions.R")
source("R/search_functions.R")

#set.seed(128)

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

lowlim=rep(-0.5,3*ls)
highlim=rep(0.5,3*ls)
#Test
search_state <- initialize_search(param_len = 3*ls,lowlim=lowlim,highlim=highlim,config = search.config(deg=1))

genstate = F
for (itt in 1:7){
  tempconfig = search.config(gen = genstate,deg=1,search_samples = ls,search_mag=0.1)
  search_state <- step_search(desired_values = ydat,search_state,ineq_w = c(),eq_w=c(),config = tempconfig)

  modout = model_function(as.vector(search_state$poly_recs[1,]),ls)
  surrogate_out <- surrogate_model(as.vector(search_state$poly_recs[1,]),search_state$polyouts,search_state$centersaves,deg=1)
  if (mean(abs((modout-surrogate_out)/modout))>0.5){
    genstate = T
  }else{
    genstate = F
  }
  matplot(xdat,search_state$current_pout,type="l",col="grey",ylim = c(0,1.2))
  points(xdat,ydat)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(search_state$poly_recs[1,]),ls),col="red")
}
#usevec = (dim(search_state$xpins)[1]-200):dim(search_state$xpins)[1]
usevec = sample(1:dim(search_state$xpins)[1],ls*2)

print(paste("Number of function evals: ",dim(search_state$xpins)[1]))

#Swap into particle swarm
swarm_state <- convert_search_to_swarm(search_state = search_state,usevec = usevec,desired_values = ydat,param_len = 3*ls,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=length(usevec),num_cluster=1,deg=1))
#swarm_state <- initialize_swarm(desired_values = ydat,param_len = 3*ls,lowlim=lowlim,highlim=highlim,ineq_w = c(),eq_w = c(),config = swarm.config(swarm_size = ls*2))

for (itt in 1:30){
  tempcontrol = swarm.control(poly_w=0.2,stoch_w = 0.01)
  swarm_state <- step_swarm(desired_values = ydat,swarm_state,ineq_w = c(),eq_w=c(),config = swarm.config(swarm_size=dim(swarm_state$x.p)[1],num_cluster=1,swarm.control = tempcontrol))
  matplot(xdat,swarm_state$pout,type="l",col="grey",ylim = c(0,1.2))
  points(xdat,ydat)
  surrogate_out <- surrogate_model(as.vector(swarm_state$poly_rec),swarm_state$polyouts,swarm_state$centersaves,deg=1)
  lines(xdat,surrogate_out,col="blue")
  lines(xdat,model_function(as.vector(swarm_state$poly_rec),ls),col="red")
}

print(paste("Number of function evals: ",dim(swarm_state$allxp)[1]))

print("-----BOBYQA-----")
bobya_function <-function(params){
  sum(abs(model_function(params)-ydat))
}
#Compare to BOBYA
result <- minqa::bobyqa(runif(3*ls,-0.5,0.5), bobya_function, lower = lowlim, upper = highlim)
lines(xdat,model_function(as.vector(result$par)),col="green")
print(result)
