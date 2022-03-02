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
ODpoly = function(powers, betas, alg = "DE", iter, swarm, pts, bound, degree = 2, crit = 'D') {
  
  # define objective function
  obj_func = obj_function_factory(powers, betas, degree, crit, bound)
  
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
  p = plot_sens(xs, out, betas, powers, degree, crit, bound)
  
  # return
  out = list(design = out, plot = p, value = c(sol$optimumValue))
}