# compute sensitivity function for a given design
# z: compute the function at z
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
sens = function(z, vars, betas, powers) {
  
  # distinguish between points and weights
  x = vars[1:pts]
  w = vars[(pts+1):(2*pts)]
  s = sum(w, na.rm = T) # need to get rid of NAs here
  if (s < 0 | s > 1) # constraint implementation
    return(-Inf)
  
  # compute eta
  # additional log terms not implemented yet
  eta = betas[1] + betas[2] * x^powers[1] + betas[3] * x^powers[2]
  etaz = betas[1] + betas[2] * z^powers[1] + betas[3] * z^powers[2]
  
  # weight functions
  sigma = exp(eta)/(1+exp(eta))^2
  sigmaz = exp(etaz)/(1+exp(etaz))^2
  
  # compute information matrix
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
  
  # compute matrix inverse and then sensitivity function
  # avoid singular matrices
  # solution from https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
  if (class(try(solve(M),silent=T))[1]!="matrix") {
    # set Minv to something
    #Minv = matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)
    y = 1
  }
  else {
    # compute sensitivity function
    Minv = solve(M)
    b = c(1, z^powers[1], z^powers[2])
    y = sigmaz * t(b) %*% Minv %*% b - 3
  }
  
  
  
  return(y)
  
  
}