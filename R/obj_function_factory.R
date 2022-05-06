# objective function factory
# generates a function for calculating the design criterion based on current solution vector
obj_function_factory = function(powers, betas, degree = 2, crit = "D", bound, p, lam) {
  
  # check input
  if (length(powers) != length(betas)-1) # make sure there is a coefficient for each power
    return(0)
  
  force(powers)
  force(betas)
  force(crit)
  force(lam)
  
  lbeta = length(betas)
  
  # define a powers vector that includes 0
  zpowers = c(0, powers)
  
  # if c-optimal design, compute gradient
  if (crit == "EDp" | crit == "Dual") {
    EDp_grad = grad_EDp(betas, powers, bound, p = p)
    if (is.na(EDp_grad$EDp))
      stop("No X value for EDp found.")
    else {
      dg = EDp_grad$grad
    }
  }
  
  # construct objective function
  # only input should be design on that point
  # vars is a list with the current x values and weights
  # switch on whether the user wants a degree 2 or 3 polynomial
  if (degree == 2) {
    
    # check correct dimensions for powers and betas
    if (length(betas) != 3 | length(zpowers) != 3) {
      print("Incorrect number of coefs or powers")
      return(0)
    }
    
    
    obj_func = function(vars) {
      
      # distinguish between points and weights
      pts = length(vars)/2
      x = vars[1:pts]
      w = vars[(pts+1):(2*pts)]
      s = sum(w, na.rm = T) # needed to fix if statement error
      if (s < 0 | s > 1) # constraint implementation
        return(-Inf)
      
      # x1 is the 2nd term in the polynomial
      x1 = H(2, x, zpowers)
      x2 = H(3, x, zpowers)
      
      # compute eta
      eta = betas[1] + betas[2] * x1 + betas[3] * x2
      
      # weight function
      sigma = exp(eta)/(1+exp(eta))^2
      
      # information matrix
      M = 0
      for (i in 1:pts) {
        
        m12 = x1[i]
        m13 = x2[i]
        m23 = x1[i]*x2[i]
        
        M_i = w[i] * sigma[i] * matrix(c(
          1, m12, m13,
          m12, x1[i]^2, m23,
          m13, m23, x2[i]^2
        ), ncol=3)
        
        M = M + M_i
        
      }
      
      # use information matrix to compute objective
      # silence warnings
      if (crit == "D") {
        obj_value = suppressWarnings(log(det(M)))
      }
      else if (crit == "A") {
        if (class(try(solve(M),silent=T))[1]!="matrix") # avoid matrix singularity
          return(Inf) # think this should also be changed to -Inf
        else {
          Minv = solve(M)
          obj_value = -sum(diag(Minv))
        }
        
      }
      else if (crit == "EDp") {
        if (class(try(solve(M),silent=T))[1]!="matrix") {
          return(-Inf)
        }
        Minv = solve(M)
        obj_value = - log(t(dg) %*% Minv %*% dg) # get as close to zero as possible
      }
      else if (crit == "Dual") {
        if (class(try(solve(M),silent=T))[1]!="matrix") {
          return(-Inf)
        }
        Minv = solve(M)
        Dval = suppressWarnings(log(det(M)))
        Cval = - suppressWarnings(log(t(dg) %*% Minv %*% dg))
        obj_value = lam * Cval + (1 - lam)/3 * Dval
      }
      else {
        obj_value = suppressWarnings(log(det(M)))
      }
      
      # deal with NAs
      if (is.na(obj_value)) 
        return(-Inf)
      else
        return(obj_value)
    }
  }
  else if (degree == 3) { # degree 3 objective function
    
    # check correct dimensions for powers and betas
    if (length(betas) != 4 | length(zpowers) != 4) {
      print("Incorrect number of coefs or powers")
      return(0)
    }
    
    obj_func = function(vars) {
      
      # distinguish between points and weights
      pts = length(vars)/2
      x = vars[1:pts]
      w = vars[(pts+1):(2*pts)]
      s = sum(w, na.rm = T) # needed to fix if statement error
      if (s < 0 | s > 1) # constraint implementation
        return(-Inf)
      
      # x1 is the 2nd term in the polynomial
      x1 = H(2, x, zpowers)
      x2 = H(3, x, zpowers)
      x3 = H(4, x, zpowers)
      
      # compute eta
      eta = betas[1] + betas[2] * x1 + betas[3] * x2 + betas[4] * x3
      
      # weight function
      sigma = exp(eta)/(1+exp(eta))^2
      
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
      
      # use information matrix to compute objective
      # using D for now
      # silence warnings
      if (crit == "D") {
        obj_value = suppressWarnings(log(det(M)))
      }
      else if (crit == "A") {
        if (class(try(solve(M),silent=T))[1]!="matrix") # avoid matrix singularity
          return(Inf)
        else {
          Minv = solve(M)
          obj_value = -sum(diag(Minv))
        }
        
      }
      else if (crit == "EDp") {
        if (class(try(solve(M),silent=T))[1]!="matrix") {
          return(-Inf)
        }
        Minv = solve(M)
        obj_value = - t(dg) %*% Minv %*% dg # get as close to zero as possible
      }
      else if (crit == "Dual") {
        if (class(try(solve(M),silent=T))[1]!="matrix") {
          return(-Inf)
        }
        Minv = solve(M)
        Dval = suppressWarnings(log(det(M)))
        Cval = - suppressWarnings(log(t(dg) %*% Minv %*% dg))
        obj_value = lam * Cval + (1 - lam)/4 * Dval
      }
      else {
        obj_value = suppressWarnings(log(det(M)))
      }
      
      # deal with NAs
      if (is.na(obj_value)) {
        if (crit == 'A')
          return(Inf)
        else
          return(-Inf)
      }
      else
        return(obj_value)
    }
  }
  else {
    print("Invalid polynomial degree")
    return(0)
  }
  
  
  return(obj_func)
}