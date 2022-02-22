# check_EDp.R
# checks all possible solutions to EDp equation and returns the feasible value and which solution was used
# if no feasible value, return NA
# p: desired ED percentile
# beta: coefficients
# powers: polynomial powers, only standard polynomials allowed
# bound: upper bound for design interval
check_EDp = function(p, beta, powers, bound) {
  
  # convert variables for algebraic sanity
  # switch on quadratic and cubic
  if (identical(powers, c(1, 2))) { # quadratic
    a = beta[3]
    b = beta[2]
    c = beta[1] - log(p/(1-p))
    
    test = b^2 - 4*a*c
    if (test < 0) # imaginary roots
      return(NA)
    else {
      # calculate roots
      x1 = (-b + sqrt(test))/(2*a)
      x2 = (-b - sqrt(test))/(2*a)
      
      # throw out roots that are out of bounds => set to Inf
      if (x1 < 0.1 | x1 > bound)
        x1 = Inf
      if (x2 < 0.1 | x2 > bound)
        x2 = Inf
      
      m = min(x1, x2) # choose lowest dose
      if (m == Inf)
        return(NA)
      else {
        if (m == x1) {
          sol = 1
        }
        else if (m == x2)
          sol = 2
        else
          sol = NA
        
        return(list(m = m, sol = sol))
        
      }
    }
  }
  else if (identical(powers, c(1, 2, 3))) { # cubic
    stop("Work in progress")
  }
  else {
    stop("Polynomial not supported")
  }
}