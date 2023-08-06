#Functions for Scedastic Surrogate Search Optimization

##############################################################################
#' Search Function configuration
#'
#' Establish parameters for surrogate search optimization
#' @param deg degree of surrogate polynomial
#' @param gen should a new search be initiated from a new location
#' @param search_mag magnitude of the local search
#' @param search_samples number of searches to perform around sample point
#' @return search configuration parameter list
#' @export
search.config <- function(deg=2,gen=F,search_mag=NULL,search_samples=NULL){
  list(deg = deg,
       gen = gen,
       search_mag = search_mag,
       search_samples = search_samples)
}

##############################################################################
#' Initialize Surrogate Search Optimization
#'
#' Function to establish initial positions of surrogate searcher
#'
#' @param param_len dimensionality of optimization problem
#' @param lowlim lower bound of input values
#' @param highlim upper bound of input values
#' @param config search configuration list
#' @return search state list
#' @export
initialize_search <- function(param_len,lowlim,highlim,config = search.config()){

  pos.x = matrix(runif(param_len,lowlim,highlim),nrow=1)
  vel =  matrix(runif(param_len,0,0),nrow=1)
  xpins = c()
  iiter = 0
  pows = NULL
  polyouts=c()
  pouts = c()

  search_state <- list(pos.x = pos.x,
                       vel = vel,
                       xpins = xpins,
                       iiter = iiter+1,
                       pows = pows,
                       polyouts = polyouts,
                       deg = config$deg,
                       highlim = highlim,
                       lowlim = lowlim,
                       pouts = c())
}

##############################################################################
#' Step through surrogate search optimization
#'
#' Function to step through surrogate search optimization by updating search states
#'
#' @param desired_vals desired outputs
#' @param search_state search state list
#' @param ineq_w inequality constraint weights
#' @param eq_w equality constraint weights
#' @param config search configuration list
#' @return search state list
#' @export
step_search <- function(desired_values,search_state,ineq_w,eq_w,config = search.config()){
  pos.x = search_state$pos.x
  vel = search_state$vel
  iiter <- search_state$iiter
  pows <- search_state$pows
  polyouts <- search_state$polyouts
  xpins <- search_state$xpins
  deg <- config$deg
  highlim <- search_state$highlim
  lowlim <- search_state$lowlim
  centersaves <- search_state$centersaves
  gen <- config$gen
  pouts <- search_state$pouts

  param_len <- dim(pos.x)[2]
  if (is.null(config$search_mag)){
    search_mag <- 0.2
  }else{
    search_mag <- config$search_mag
  }

  if (is.null(config$search_samples)){
    search_samples <- param_len*2
  }else{
    search_samples <- config$search_samples
  }

  vel[(pos.x+vel)<lowlim] = -vel[(pos.x+vel)<lowlim]/100
  vel[(pos.x+vel)>highlim] = -vel[(pos.x+vel)>highlim]/100
  pos.x_old <- pos.x
  pos.x <- pos.x+vel

  if (gen){
    pos.x = rbind(pos.x,runif(param_len,lowlim,highlim))
  }

  poly_recs = c()
  current_pout = c()
  for (ip in sample(1:dim(pos.x)[1])){
    centersaves = rbind(centersaves,pos.x[ip,])

    search_low = pmax(lowlim,pos.x[ip,]-(highlim-lowlim)*search_mag)
    search_high = pmin(highlim,pos.x[ip,]+(highlim-lowlim)*search_mag)

    xpin = rbind(sample_lhs_around_center(pos.x[ip,],
                             n_samples = search_samples,
                             search_low-pos.x[ip,],
                             search_high-pos.x[ip,]),
                     pos.x[ip,])
    pout <- apply(xpin,1,modf)
    current_pout <- cbind(current_pout,pout[,dim(pout)[2]])
    pouts <- cbind(pouts,pout)

    outgs<-apply(pout,2,function(aa)fitness_fun(aa,desired_values,ineq_w,eq_w))

    permat = t(t(xpin)-pos.x[ip,])
    X_poly <- make_poly(permat, degree=deg)

    if (iiter==1){
      if (deg==1){
        pows = diag(param_len)
      }else{
        pows = apply(do.call(rbind,strsplit(names(X_poly[1,]),"[.]")),2,as.numeric)
      }
    }

    temp_polys=t(apply(pout,1,function(oo)fit_poly(X_poly,oo)))
    polyouts=abind::abind(polyouts,temp_polys,along=3)

    xpins = rbind(xpins,xpin)

    control <- list(NP = min(param_len,100),
                    itermax = 30,
                    F = 0.8,
                    CR = 0.9,
                    trace=FALSE
    )

    result <- suppressWarnings(DEoptim::DEoptim(poly_obj_fn,
                                                lower=search_low,
                                                upper=search_high,
                                                control,
                                                polys = polyouts,
                                                desired_values = desired_values,
                                                centers=centersaves,
                                                deg=deg))
    poly_rec = array(result$optim$bestmem)
    poly_recs = rbind(poly_recs,poly_rec)
  }

  vel = 1*(poly_recs - pos.x)

  search_state <- list(pos.x = pos.x,
                       vel = vel,
                       xpins = xpins,
                       iiter = iiter+1,
                       pows = pows,
                       polyouts = polyouts,
                       deg = search_state$deg,
                       centersaves=centersaves,
                       highlim=highlim,
                       lowlim=lowlim,
                       poly_recs = poly_recs,
                       pouts = pouts,
                       current_pout=current_pout)
}

##############################################################################
#' Perform Latin Hypercube Sampling around a Center Point
#'
#' Function to do Latin Hypercube Sampling around a central location based on
#' allowable ranges
#'
#' @param center center point location
#' @param n_samples number of Latin Hypercube samples to perform
#' @param range_min lower bound of sampling
#' @param range_max upper bound of sampling
#' @return latin hypercube samples
#' @export
sample_lhs_around_center <- function(center, n_samples, range_min, range_max) {
  n_dim <- length(center)

  lhs_samples <- lhs::randomLHS(n = n_samples, k = n_dim)
  #lhs_samples <- lhs::optimumLHS(n_samples,n_dim,10,0.05)

  scaled_samples <- sweep(lhs_samples, 2, range_max - range_min, FUN = "*")
  shifted_samples <- sweep(scaled_samples, 2, range_min, FUN = "+")
  centered_samples <- sweep(shifted_samples, 2, center, FUN = "+")

  return(centered_samples)
}

##############################################################################
#' Convert Search Results into Initial Scedastic Surrogate Swarm Position
#'
#' Function to convert the results of a surrogate search into the initial locations
#' of a Scedastic Surrogate Swarm optimization, preserving the learned behaviour
#' by passing saved surrogates
#'
#' @param search_state search state list
#' @param usevec which of the search instances to use
#' @param desired_vals desired outputs
#' @param param_len dimensionality of optimization problem
#' @param lowlim lower bound of input values
#' @param highlim upper bound of input values
#' @param ineq_w inequality constraint weights
#' @param eq_w equality contstraint weights
#' @param config swarm configuration list
#' @return swarm state list
#' @export
convert_search_to_swarm <- function(search_state,usevec,desired_values,param_len,lowlim,highlim,ineq_w,eq_w,config = swarm.config()){
  numineq = length(ineq_w)
  numeq = length(eq_w)

  x.p = search_state$xpins[usevec,]
  pout = search_state$pout[,usevec]

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

  search_state <- list(x.p=x.p,
                      vel=vel,
                      best_p=best_p,
                      bestgs=bestgs,
                      errbest=errbest,
                      locality=locality,
                      phet=phet,
                      pout=pout,
                      allxp=search_state$xpins,
                      allouts=search_state$pout,
                      centersaves=search_state$centersaves,
                      polyouts=search_state$polyouts,
                      clustrack=c(),
                      iiter = 1,
                      deg = config$deg,
                      highlim = highlim,
                      lowlim = lowlim,
                      pows = NULL)
}
