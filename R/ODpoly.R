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
ODpoly = function(powers, betas, alg = "DE", iter, swarm, pts, bound) {
  
  # define objective function
  obj_func = obj_function_factory(powers, betas)
  
  # set up for metaheuristics
  numVar = 2*pts # each point has one weight and last weight is sum of the others
  d = c(rep(c(0, bound), pts), rep(c(0,1), pts))
  rangeVar = matrix(d, nrow = 2)
  
  # algorithm params
  control = list(numPopulation = swarm, maxIter = iter)
  
  # find optimal design
  sol = metaOpt(obj_func, optimType = "MAX", algorithm = alg, 
                                  numVar, rangeVar, control)
  
  # save output
  result = sol$result
  
  # plot sensitivity function
  xs = seq(0, bound, 0.01)
  p = plot_sens(xs, sol$result, betas, powers)
  
  # return
  out = list(design = result, plot = p)
}