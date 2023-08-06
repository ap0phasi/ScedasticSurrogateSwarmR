#Functions for Scedastic Surrogate Swarm Optimization

##############################################################################
#' Construct Polynomial Instance
#'
#' This function creates a polynomial matrix of the specified degree provided
#' input data
#'
#' @param dat input data
#' @param degree polynomial degree
#' @return polynomial input array
#' @export
make_poly <- function(dat,degree){
  if (degree==1){
    return(dat)
  }else{
    return(poly(dat, degree=degree, raw = TRUE,simple=TRUE))
  }
}

##############################################################################
#' Customization of error function
#'
#' Allows for different error functions to be used
#'
#' @param actual known values
#' @param predicted desired values
#' @param type what error function should be used?
##' \itemize{
##'  \item{"abs"}{ Sum of absolute errors}
##'  \item{"rme"}{ Root mean squared error}
##' }
#' @return error value
error_function <- function(actual,predicted,type="abs"){
  if (type=="rmse") {
    error <- sqrt(mean((actual - predicted)^2))
  }else if (type=="abs"){
    error <- sum(abs(actual - predicted))
  } else {
    error <- mean(abs(actual - predicted))
  }
  return(error)
}

##############################################################################
#' Fit Polynomial Instance
#'
#' Function to fit polynomial to observed data
#'
#' @param X_poly Polynomial input array
#' @param y Response values to fit polynomial to
#' @param crossval Should cross validation be used in fitting?
#' @return polynomial coefficients
#' @export
fit_poly <- function(X_poly, y,crossval = F) {
  if (length(y)<100){
    poly_fit <- lm(y ~ X_poly)
    poly_coef <- coef(poly_fit)
    poly_coef[is.na(poly_coef)]=0
  }else{
    if (crossval){
      cvfit <- glmnet::cv.glmnet(X_poly, y,maxit=10^3)
      poly_coef <- as.numeric(coef(cvfit, s = "lambda.min"))
    }else{
      poly_fit <- glmnet::glmnet(X_poly, y,maxit=10^3)
      poly_coef <- as.numeric(coef(poly_fit, s = 0.001))
    }
  }
  return(poly_coef)
}

##############################################################################
#' Generalized Normalization Function
#'
#' Function to perform robust normalization in case the sum of x is 0.
#'
#' @param x unnormalized array
#' @return normalized array
normalize <- function(x){
  if (!sum(x)==0){
    return(x/sum(x))
  }else{
    return(x)
  }
}

##############################################################################
#' Fitness Function
#'
#' Evaluate accuracy of solution
#'
#' @param mod_vals outputs of model function with current inputs
#' @param desired_vals desired outputs
#' @param ineq_w weights of inequality constraints
#' @param eq_w weights of equality constraints
#' @return fitness score
#' @export
fitness_fun <- function(mod_vals,desired_vals,ineq_w,eq_w){
  numineq = length(ineq_w)
  numeq = length(eq_w)
  l_mod_val <- length(desired_vals)-numineq-numeq
  mod_score <- error_function(mod_vals[1:l_mod_val],desired_vals[1:l_mod_val])
  if (numeq>0){
    eq_vals <- mod_vals[(l_mod_val+1):(numeq+l_mod_val)]
    eq_score <- sum(abs(eq_vals))
  }else{
    eq_score=c()
  }
  if (numineq>0){
    ineq_vals <- mod_vals[(l_mod_val+numeq+1):(numineq+numeq+l_mod_val)]
    ineq_score <- sum(pmax(ineq_vals*0.1,ineq_vals))
  }else{
    ineq_score <- c()
  }
  return(sum(c(mod_score,eq_score*eq_w,ineq_score*ineq_w)))
}

##############################################################################
#' Evaluate Surrogate Model
#'
#' Evaluate the constructed surrogate model of the true function based on proximity
#' to centers and respective polynomials fit to response
#' @param xsel selected input values to evaluate surrogate for
#' @param polys all polynomials fit to explored responses
#' @param centers polynomial centers
#' @param deg degree of polynomial (recommended between 0 and 1)
#' @param type should we only use the polynomial evaluated with the center closest to xsel?
##' \itemize{
##'  \item{"closest"}{ Only use the polynomial whose center is closest to xsel}
##'  \item{"all"}{ Evaluate all polynomials and do weighting average based on distance of respective centers to xsel}
##' }
#' @param form should distance be evaluated with a simple absolute distance?
##' \itemize{
##'  \item{"matrix"}{ Evaluate distance as sum of absolute differences}
##'  \item{"errorfun"}{ Evaluate distance with error function}
##' }
#' @return predicted model outputs based on surrogate
#' @export
surrogate_model <- function(xsel, polys, centers, deg,type = "closest",form = "matrix"){
  if (type=="closest"){
    if (form=="matrix"){
      closest_center <- which.min(colSums(abs(xsel-t(centers))))[1]
    }else{
      closest_center <- which.min(apply(centers,1,function(cc)error_function(cc,xsel,type = "abs")))[1]
    }
    polysel = cbind(1,make_poly(t(xsel-t(centers)),degree=deg))
    modeled_values <- as.vector(polys[,,closest_center]%*%array(polysel[closest_center,]))
  }else{
    if (form=="matrix"){
      closest_vals <- order(colSums(abs(xsel-t(centers))))[1:min(2,dim(centers)[1])]
      polysel = cbind(1,make_poly(t(xsel-t(centers)),degree=deg))
      modeled_array = c()
      for (icc in closest_vals){
        modeled_array <- rbind(modeled_array,as.vector(polys[,,icc]%*%array(polysel[icc,])))
      }
      modeled_values <- t(modeled_array)%*%normalize((length(closest_vals):1)^50)
    }
  }
  return(modeled_values)
}

##############################################################################
#' Objective Function for Suboptimization Problem
#'
#' Construct objective function for the suboptimization
#' @param xsel selected input values to evaluate surrogate for
#' @param polys all polynomials fit to explored responses
#' @param desired_vals desired outputs
#' @param centers polynomial centers
#' @param deg degree of polynomial (recommended between 0 and 1)
#' @return objective function score
#' @export
poly_obj_fn <- function(xsel, polys, desired_values,centers,deg) {
  modeled_values <- surrogate_model(xsel,polys,centers,deg)
  obj_value <- error_function(desired_values,modeled_values)
  return(obj_value)
}

##############################################################################
#' Heteroscedastic Loss Function
#'
#' Function to determine heteroscedastic loss
#' @param goal desired outputs
#' @param mean modeled means
#' @param sd modeled standard deviations
#' @return heteroscedastic loss
#' @export
hetloss <- function(goal,mean,sd){
  return(1/(2*sd^2)*abs(goal-mean)^2+1/2*log(sd^2))
}


##############################################################################
#' Cluster Particles in Swarm
#'
#' Group particles in swarm into clusters
#' @param x.p positions of particles in swarm
#' @param num_cluster how many clusters should the particles be grouped into?
#' @param min_points what is the minimum number of partciles that should be in a cluster?
#' @return cluster assignment for each particle in swarm
#' @return cluster centers
#' @export
cluster_swarm <- function(x.p,num_cluster,min_points){
  cluster_res = kmeans(x.p,num_cluster)
  small_clusters <- names(table(cluster_res$cluster)[table(cluster_res$cluster) < min_points])
  for (i in small_clusters) {
    distances <- apply(matrix(x.p[cluster_res$cluster == i, ],ncol=dim(x.p)[2]), 1, function(x) sqrt(rowSums(t(x-t(cluster_res$centers))^2)))
    new_cluster <- apply(distances,2,function(d)c(1:length(d))[-as.numeric(i)][which.min(d[-as.numeric(i)])])
    cluster_res$cluster[cluster_res$cluster == i] <- new_cluster
  }
  return(cluster_res)
}

##############################################################################
#' Scedastic Surrogate Swarm Control
#'
#' Establish weights and importances for Scedastic Surrogate Swarm for determining new
#' particle velocities.
#' @param stoch_w how much influence does scedastic factor have?
#' @param vel_w how much influence does the previous velocity have?
#' @param r_p_w how much influence does the particle's historical best position have?
#' @param r_g_w how much influence does the global current best position have?
#' @param r_go_w how much influence does the global historical best position have?
#' @param poly_w how much influence does the suboptimized surrogate function have?
#' @param sens_overall do we use the overall sensitivity when backpropagating the stochastic factor?
#' @return control list
#' @export
swarm.control <- function(stoch_w = 0.1,
                      vel_w = 0.5,
                      r_p_w = 0.3,
                      r_g_w = 0.6,
                      r_go_w = 0.6,
                      poly_w = 0.2,
                      sens_overall = F){
  list(stoch_w = stoch_w,
       vel_w = vel_w,
       r_p_w = r_p_w,
       r_g_w = r_g_w,
       r_go_w = r_go_w,
       poly_w = poly_w,
       sens_overall = sens_overall)
}

##############################################################################
#' Scedastic Surrogate Swarm Configuration
#'
#' Establish parameters for evaluating Scedastic Surrogate Swarm
#' @param swarm_size number of particles in swarm
#' @param num_cluster number of clusters for grouping particles
#' @param min_points minimum number of particles per cluster
#' @param deg degree of surrogate polynomial
#' @param localityfac what level of locality to use when determining global values (1 would be full global topology)
#' @param opt_from_best when performing suboptimization should we start from the best known position?
#' @param swarm.control swarm control weight list
#' @return swarm configuration parameter list
#' @export
swarm.config <- function(swarm_size = 100,
                         num_cluster = 1,
                         min_points = 3,
                         deg = 1,
                         localityfac=0.6,
                         opt_from_best = F,
                         swarm.control = list(stoch_w = 0.1,
                                              vel_w = 0.5,
                                              r_p_w = 0.3,
                                              r_g_w = 0.6,
                                              r_go_w = 0.6,
                                              poly_w = 0.2,
                                              sens_overall = F)){
  list(swarm_size = swarm_size,num_cluster = num_cluster,min_points = min_points,deg = deg,localityfac=localityfac,opt_from_best=opt_from_best,swarm.control=swarm.control)
}

##############################################################################
#' Initialize Scedastic Surrogate Swarm Optimization
#'
#' Function to establish initial positions and velocities of Scedastic Surrogate Swarm
#'
#' @param desired_vals desired outputs
#' @param param_len dimensionality of optimization problem
#' @param lowlim lower bound of input values
#' @param highlim upper bound of input values
#' @param ineq_w inequality constraint weights
#' @param eq_w equality contstraint weights
#' @param config swarm configuration list
#' @return swarm state list
#' @export
initialize_swarm <-function(desired_values,param_len,lowlim,highlim,ineq_w,eq_w,config = swarm.config()){
  numineq = length(ineq_w)
  numeq = length(eq_w)
  x.p <- matrix(runif(param_len*config$swarm_size,
                   rep(lowlim,config$swarm_size),
                   rep(highlim,config$swarm_size)),nrow=config$swarm_size)
  pout <- apply(x.p,1,modf)
  pmean <- apply(pout,1,mean)
  psd <- apply(pout,1,sd)
  phet <- hetloss(desired_values,pmean,psd)
  perr <- abs(pout-desired_values)

  errbest <- perr
  outgs<-apply(pout,2,function(aa)fitness_fun(aa,desired_values,ineq_w,eq_w))
  bestgs<-outgs
  best_p<-x.p

  pbest<-replicate(length(desired_values),x.p)
  gbest<-x.p[apply(perr,1,function(x) which.min(x)),]
  gerrbest<-apply(perr,1,min)

  vel <- matrix(runif(param_len*config$swarm_size,-0.1,0.1),ncol=param_len)

  locality <- round(config$localityfac*config$swarm_size)

  swarm_state <- list(x.p=x.p,
                      vel=vel,
                      best_p=best_p,
                      bestgs=bestgs,
                      errbest=errbest,
                      locality=locality,
                      phet=phet,
                      pout=pout,
                      allxp=c(),
                      allouts=c(),
                      centersaves=c(),
                      polyouts=c(),
                      clustrack=c(),
                      iiter = 1,
                      deg = config$deg,
                      highlim = highlim,
                      lowlim = lowlim,
                      pows = NULL)
}

##############################################################################
#' Step through Scedastic Surrogate Swarm optimization
#'
#' Function to step through Scedastic Surrogate Swarm optimization by updating swarm states
#'
#' @param desired_vals desired outputs
#' @param swarm_state swarm state list
#' @param ineq_w inequality constraint weights
#' @param eq_w equality constraint weights
#' @param config swarm configuration list
#' @return swarm state list
#' @export
step_swarm <- function(desired_values,swarm_state,ineq_w,eq_w,config = swarm.config()){
  x.p <- swarm_state$x.p
  vel <- swarm_state$vel
  best_p <- swarm_state$best_p
  bestgs <- swarm_state$bestgs
  errbest <- swarm_state$errbest
  locality <- swarm_state$locality
  phet <- swarm_state$phet
  allxp <- swarm_state$allxp
  allouts <- swarm_state$allouts
  centersaves <- swarm_state$centersaves
  polyouts <- swarm_state$polyouts
  clustrack <- swarm_state$clustrack
  iiter <- swarm_state$iiter
  deg <- swarm_state$deg
  highlim <- swarm_state$highlim
  lowlim <- swarm_state$lowlim
  pows <- swarm_state$pows

  numineq = length(ineq_w)
  numeq = length(eq_w)


  swarm_size <- dim(x.p)[1]
  param_len <- dim(x.p)[2]

  locality <- round(config$localityfac*swarm_size)

  hetold <- phet
  vel[(x.p+vel)<lowlim] = -vel[(x.p+vel)<lowlim]/100
  vel[(x.p+vel)>highlim] = -vel[(x.p+vel)>highlim]/100
  x.p <- x.p+vel

  pout <- apply(x.p,1,modf)
  pmean <- apply(pout,1,mean)
  psd <- apply(pout,1,sd)
  phet <- hetloss(desired_values,pmean,psd)
  perr <- abs(pout-desired_values)
  allxp <- rbind(allxp,x.p)
  allouts <- cbind(allouts,pout)

  outgs<-apply(pout,2,function(aa)fitness_fun(aa,desired_values,ineq_w,eq_w))
  cmat<-t(apply(as.matrix(dist(x.p,method="manhattan")),1,function(x)order(x)[-1][1:locality]))
  best_g_mat <- t(apply(cmat,1,function(a) x.p[a,][which.min(outgs[a]),]))
  best_go_mat <- t(apply(cmat,1,function(a) best_p[a,][which.min(bestgs[a]),]))
  new.ind=which(outgs<bestgs)
  best_p[new.ind,]<-x.p[new.ind,]
  bestgs[new.ind]<-outgs[new.ind]

  r_p=matrix(runif(swarm_size*param_len,0,1),nrow=swarm_size,ncol=param_len)
  r_g=matrix(runif(swarm_size*param_len,0,1),nrow=swarm_size,ncol=param_len)
  r_go=matrix(runif(swarm_size*param_len,0,1),nrow=swarm_size,ncol=param_len)
  r_poly=matrix(runif(swarm_size*param_len,0,1),nrow=swarm_size,ncol=param_len)

  cluster_res <- cluster_swarm(x.p,config$num_cluster,config$min_points)


  rback = c()
  for (icc in 1:dim(cluster_res$centers)[1]){
    clusel=x.p[which(cluster_res$cluster==icc),]
    xaps = which.min(rowSums(abs(t(t(allxp)-cluster_res$centers[icc,]))))
    centersaves=rbind(centersaves,allxp[xaps,])
    xpin = rbind(x.p[cluster_res$cluster==icc,],0)
    permat = t(t(xpin)-cluster_res$centers[icc,])

    X_poly <- make_poly(permat, degree=deg)

    if (iiter==1){
      if (deg==1){
        pows = diag(param_len)
      }else{
        pows = apply(do.call(rbind,strsplit(names(X_poly[1,]),"[.]")),2,as.numeric)
      }
    }

    temp_polys=t(apply(cbind(pout[,cluster_res$cluster==icc],allouts[,xaps]),1,function(oo)fit_poly(X_poly,oo)))
    polyouts=abind::abind(polyouts,temp_polys,along=3)

    clus_size = length(which(cluster_res$cluster==icc))
    if (clus_size>0&!config$swarm.control$sens_overall){
      sens = t(apply(temp_polys,1,function(aa)t(t(pows)%*%abs(t(t(aa[-1]))))))
      backprop = array(MASS::ginv(sens)%*%pmax(0,phet))
      rmag = abs(backprop/clus_size)
      rback = rbind(rback,t(qnorm(t(lhs::randomLHS(clus_size,param_len)),0,rmag)))
      #rback = rbind(rback,t(qnorm(t(qunif(lhs::randomLHS(clus_size,param_len),0.1,0.9)),0,rmag)))
    }
  }

  control <- list(NP = min(param_len,100),
                  itermax = 30,
                  F = 0.8,
                  CR = 0.9,
                  trace=FALSE
  )

  if (config$opt_from_best){
    start=best_p[sample(1:swarm_size,1),]
  }else{
    start=x.p[sample(1:swarm_size,1),]
  }

  result <- suppressWarnings(DEoptim::DEoptim(poly_obj_fn,
                                              lower=start-0.3*abs(start),
                                              upper=start+0.3*abs(start),
                                              control,
                                              polys = polyouts,
                                              desired_values = desired_values,
                                              centers=centersaves,
                                              deg=deg))

  poly_rec = array(result$optim$bestmem)

  if (config$swarm.control$sens_overall){
    sens_weights=normalize((1:dim(polyouts)[3])^100)
    sens = t(apply(polyouts[,,1],1,function(aa)t(t(pows)%*%abs(t(t(aa[-1]))))))*sens_weights[1]
    if (dim(polyouts)[3]>1){
      for (iss in 2:dim(polyouts)[3]){
        sens = sens+t(apply(polyouts[,,iss],1,function(aa)t(t(pows)%*%abs(t(t(aa[-1]))))))*sens_weights[iss]
      }
    }

    backprop = array(MASS::ginv(sens)%*%pmax(0,phet))
    rback = t(replicate(swarm_size,rnorm(param_len,0,abs(backprop/swarm_size))))
    #rmag = abs(backprop/swarm_size)
    #rback = t(qnorm(t(lhs::randomLHS(swarm_size,param_len)),0,rmag))
    #rback = t(qnorm(t(qunif(lhs::randomLHS(swarm_size,param_len),0.1,0.9)),0,rmag))
  }


  vel=rback*config$swarm.control$stoch_w+
    vel*config$swarm.control$vel_w+
    r_p*config$swarm.control$r_p_w*(best_p-x.p)+
    r_g*config$swarm.control$r_g_w*(best_g_mat-x.p)+
    r_go*config$swarm.control$r_go_w*(best_go_mat-x.p)+
    r_poly*config$swarm.control$poly_w*(t(matrix(poly_rec))[rep(1,swarm_size),]-x.p)

  swarm_state <- list(x.p=x.p,
                      vel=vel,
                      best_p=best_p,
                      bestgs=bestgs,
                      errbest=errbest,
                      locality=locality,
                      phet=phet,
                      pout=pout,
                      allxp=allxp,
                      allouts=allouts,
                      centersaves=centersaves,
                      polyouts=polyouts,
                      clustrack=clustrack,
                      iiter = iiter + 1,
                      deg = swarm_state$deg,
                      highlim = highlim,
                      lowlim = lowlim,
                      pows = pows,
                      poly_rec = poly_rec
  )
}

##############################################################################
#' Model Function
#'
#' Concatenate outputs of user defined model function with user defined inequality
#' and equality constraints
#'
#' @param params model input array
#' @return model output array with constraint values
#' @export
modf <- function(params){
  modval <- model_function(params)
  modout <- c(modval,eq(params),ineq(params))
  return(modout)
}
