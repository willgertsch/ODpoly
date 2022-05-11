# ODpoly
# finds optimal design for binary outcome with polynomial linear predictor
# this is a wapper that makes it easy for the user
# inputs:
# powers: 2 vector of powers for degree 2 fractional poly
# betas: 3 vector of coefficients
# alg: a metaheuristic from metaheuristicOpt
# iter: number of iterations
# swarm: size of swarm
# pts: number of design points
# bound: upper bound for design interval
# crit: design criterion to use
# p: percentile for EDp designs
# lam: weight parameter for dual objective designs
library(metaheuristicOpt)
ODpoly = function(powers, betas, alg = "DE", iter, swarm, pts, bound,
                  crit = 'D', p = 0.5, lam = 0) {
  
  degree = length(powers)
  
  # define objective function
  obj_func = obj_function_factory(powers, betas, crit, bound, p, lam)
  
  # separate flow for using RccpDE package
  
  # set up for metaheuristics
  numVar = 2*pts # each point has one weight and last weight is sum of the others
  d = c(rep(c(0.1, bound), pts), rep(c(0,1), pts)) # lower bound tries to be close to 0
  rangeVar = matrix(d, nrow = 2)
  
  # algorithm params
  control = list(numPopulation = swarm, maxIter = iter)
  
  # find optimal design
  sol = metaOpt(obj_func, optimType = "MAX", algorithm = alg, 
                numVar, rangeVar, control)
  
  
  # save output
  result = sol$result
  
  # combine design points
  raw = result
  l = length(raw)
  
  # purge points with zero weight
  # out to 3 decimal points to match output
  if (sum(round(raw[(l/2 + 1):l], 3)==0) > 0) {
    x_indices = which(round(raw[(l/2 + 1):l], 3)==0)
    w_indices = x_indices + l/2
    raw = raw[,-c(x_indices, w_indices)]
    l = length(raw)
  }
  
  # combine weights of identical points
  # sort as well
  xs = raw[1:(l/2)]
  ws = raw[(l/2+1):l]
  if (length(unique(xs)) != length(xs)) {
    dups = xs[duplicated(xs)] # keep track of dups
    for (d in dups) { # iterate through and combine weights
      indices = xs == d
      new_w = sum(ws[indices]) # add w's for a specific duplicate
      ws[indices] = new_w # update w's; will drop later
    }
    xs = unique(xs) # remove duplicates
    ws = unique(ws)
    
    raw = c(xs, ws)
    l = length(raw)
  }
  
  # labeling
  # probably a better way to do this
  # labs is a function => call it labbs
  labbs = character(l)
  for (i in 1:(l/2)) {
    labbs[i] = paste("x", toString(i), sep = "")
  }
  for (i in (l/2 + 1):l) {
    labbs[i] = paste("w", toString(i-l/2), sep = "")
  }
  
  # magic
  raw = c(raw)
  
  # sort by x's
  raw_x = raw[1:(l/2)]
  raw_w = raw[(l/2 + 1):l]
  r = rank(raw_x)
  raw_x = raw_x[r]
  raw_w = raw_w[r]
  raw = c(raw_x, raw_w)
  
  names(raw) = labbs
  
  out = raw

  # plot sensitivity function
  step = bound/1000
  xs = seq(0.1, bound, step)
  p = plot_sens(xs, out, betas, powers, crit, bound, p, lam)
  
  # return
  out = list(design = out, plot = p, value = c(sol$optimumValue))
}