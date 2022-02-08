# ODpoly
# finds optimal design for binary outcome with polynomial linear predictor
# this is a wapper that makes it easy for the user
# hopefully the user won't have to mess around with algorithm choices
# inputs:
# powers: 2 vector of powers for degree 2 fractional poly
# betas: 3 vector of coefficients
# alg: a metaheuristic from metaheuristicOpt
# iter: number of iterations
# swarm: size of swarm
# pts: number of design points
# bound: upper bound for design interval
library(metaheuristicOpt)
library(DEoptim)
ODpoly = function(powers, betas, alg = "DE", iter, swarm, pts, bound, degree = 2) {
  
  # define objective function
  obj_func = obj_function_factory(powers, betas, degree)
  
  # separate flow for using RccpDE package
  #cat(file=stderr(), length(alg))
  if (alg == "FDE") {
    # have to flip objective function
    neg_obj_func = function(x) {
      return(-obj_func(x))
    }
    
    # set bounds
    lower = c(rep(0.1, pts), rep(0, pts))
    upper = c(rep(bound, pts), rep(1, pts))
    
    # initial values
    # start as equally weighted uniform design
    # init = matrix(rep(c(seq(.1, bound, length.out=pts),
    #                     rep(1/pts, pts)), swarm), nrow=swarm, byrow=T)
    
    # equally weighted random design
    # init = matrix(rep(c(runif(pts, 0.1, bound),
    #                     rep(1/pts, pts)), swarm), nrow=swarm, byrow=T)
    
    # set algorithm params
    control = DEoptim.control(NP = swarm, itermax = iter, trace = F)
    
    # run faster implementation of DE
    out = DEoptim(neg_obj_func, lower, upper, control)
    sol = list(result = out$optim$bestmem, optimumValue = -out$optim$bestval)
  }
  else {
    # set up for metaheuristics
    numVar = 2*pts # each point has one weight and last weight is sum of the others
    d = c(rep(c(0.1, bound), pts), rep(c(0,1), pts)) # lower bound tries to be close to 0
    rangeVar = matrix(d, nrow = 2)
    
    # algorithm params
    control = list(numPopulation = swarm, maxIter = iter)
    
    # find optimal design
    sol = metaOpt(obj_func, optimType = "MAX", algorithm = alg, 
                  numVar, rangeVar, control)
  }
  
  # save output
  result = sol$result
  
  
  # plot sensitivity function
  step = bound/1000
  xs = seq(0.1, bound, step)
  p = plot_sens(xs, sol$result, betas, powers, degree)
  
  # return
  out = list(design = result, plot = p, value = c(sol$optimumValue))
}