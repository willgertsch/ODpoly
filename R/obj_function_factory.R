# objective function factory
# generates a function for calculating the design criterion based on current solution vector
obj_function_factory = function(powers, betas) {
  
  # check input
  if (length(powers) != length(betas)-1) # make sure there is a coefficient for each power
    return(0)
  
  force(powers)
  force(betas)
  
  lbeta = length(betas)
  
  # construct objective function
  # only input should be design on that point
  # vars is a list with the current x values and weights
  obj_func = function(vars) {
    
    # distinguish between points and weights
    x = vars[1:pts]
    w = vars[(pts+1):(2*pts)]
    s = sum(w)
    if (s < 0 | s > 1) # constraint implementation
      return(-Inf)
    
    # check if x > 0 for fractional poly
    
    
    # compute eta
    # additional log terms not implemented yet
    eta = betas[1] + betas[2] * x^powers[1] + betas[3] * x^powers[2]
    
    
    # weight function
    sigma = exp(eta)/(1+exp(eta))^2
    
    # information matrix
    # currently quadratic
    M = 0
    for (i in 1:pts) {
      
      # will need to update to use correct x functions
      m12 = x[i]^powers[1]
      m13 = x[i]^powers[2]
      m23 = x[i]^(powers[1]+powers[2])
      
      M_i = w[i] * sigma[i] * matrix(c(
        1, m12, m13,
        m12, x[i]^(2*powers[1]), m23,
        m13, m23, x[i]^(2*powers[2])
      ), ncol=3)
      
      M = M + M_i
      
    }
    
    
    # use information matrix to compute objective
    # using D for now
    # silence warnings
    
    obj_value = suppressWarnings(log(det(M)))
    
    
    
    # deal with NAs
    if (is.na(obj_value)) {
      return(-Inf)
    }
    else
      return(obj_value)
  }
  
  return(obj_func)
}