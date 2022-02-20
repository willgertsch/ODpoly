# compute sensitivity function for a given design
# z: compute the function at z
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
sens = function(z, vars, betas, powers, degree = 2, crit = "D") {
  
  # switch for different degree polynomials
  if (degree == 2) {
    # distinguish between points and weights
    pts = length(vars)/2
    x = vars[1:pts]
    w = vars[(pts+1):(2*pts)]
    s = sum(w, na.rm = T) # need to get rid of NAs here
    if (s < 0 | s > 1) # constraint implementation
      return(-Inf)
    
    # define a powers vector that includes 0
    zpowers = c(0, powers)
    
    x1 = H(2, x, zpowers)
    x2 = H(3, x, zpowers)
    z1 = H(2, z, zpowers)
    z2 = H(3, z, zpowers)
    
    # compute eta
    eta = betas[1] + betas[2] * x1 + betas[3] * x2
    etaz = betas[1] + betas[2] * z1 + betas[3] * z2
    
    # weight functions
    sigma = exp(eta)/(1+exp(eta))^2
    sigmaz = exp(etaz)/(1+exp(etaz))^2
    
    # compute information matrix
    M = 0
    for (i in 1:pts) {
      
      m12 = x1[i]
      m13 = x2[i]
      m23 = x1[i] * x2[i]
      
      M_i = w[i] * sigma[i] * matrix(c(
        1, m12, m13,
        m12, x1[i]^2, m23,
        m13, m23, x2[i]^2
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
      Minv = solve(M)
      b = c(1, z1, z2)
      # compute sensitivity function
      if (crit == "D") {
        y = sigmaz * t(b) %*% Minv %*% b - 3
      }
      else if (crit == "A") {
        Minv2 = Minv %*% Minv
        y = sigmaz * t(b) %*% Minv2 %*% b - sum(diag(Minv))
      }
      else { # default to D
        y = sigmaz * t(b) %*% Minv %*% b - 3
      }
      
    }
  }
  else if (degree == 3) {
    # distinguish between points and weights
    pts = length(vars)/2
    x = vars[1:pts]
    w = vars[(pts+1):(2*pts)]
    s = sum(w, na.rm = T) # need to get rid of NAs here
    if (s < 0 | s > 1) # constraint implementation
      return(-Inf)
    
    # define a powers vector that includes 0
    zpowers = c(0, powers)
    
    x1 = H(2, x, zpowers)
    x2 = H(3, x, zpowers)
    x3 = H(4, x, zpowers)
    z1 = H(2, z, zpowers)
    z2 = H(3, z, zpowers)
    z3 = H(4, z, zpowers)
    
    # compute eta
    eta = betas[1] + betas[2] * x1 + betas[3] * x2 + betas[4] * x3
    etaz = betas[1] + betas[2] * z1 + betas[3] * z2 + betas[4] * z3
    
    # weight functions
    sigma = exp(eta)/(1+exp(eta))^2
    sigmaz = exp(etaz)/(1+exp(etaz))^2
    
    # information matrix
    M = 0
    for (i in 1:pts) {
      
      m12 = x1[i]
      m13 = x2[i]
      m14 = x3[i]
      m23 = x1[i]*x2[i]
      m24 = x1[i]*x3[i]
      m34 = x2[i]*x3[i]
      
      M_i = w[i] * sigma[i] * matrix(c(
        1, m12, m13, m14,
        m12, x1[i]^2, m23, m24,
        m13, m23, x2[i]^2, m34,
        m14, m24, m34, x3[i]^2
      ), ncol=4)
      
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
      b = c(1, z1, z2, z3)
      if (crit == "D") {
        y = sigmaz * t(b) %*% Minv %*% b - 4
      }
      else if (crit == "A") {
        Minv2 = Minv %*% Minv
        y = sigmaz * t(b) %*% Minv2 %*% b - sum(diag(Minv))
      }
      else { # default to D
        y = sigmaz * t(b) %*% Minv %*% b - 4
      }
      
    }
  }
  else {
    print("Polynomial degree not supported")
    return(-1)
  }
  
  # return computed sensitivity function
  return(y)
  
}