# EDp_grad.R
# compute the gradient for EDp  to find c-optimal designs
# returns p x 1 gradient evaluated at current model parameters
# support quadratic only for now
# model: values are "quadratic", "cubic", "fp"
# solution: which solution to take gradient of
EDp_grad = function(beta, powers, solution) {
  
  if (identical(powers, c(1, 2))) { # quadratic
    # translate parameters
    c = beta[1]
    b = beta[2]
    a = beta[3]
    
    if (solution == 1) {
      dbeta1 = -1 / sqrt(b^2 - 4 * a * c)
      dbeta2 = (b/sqrt(b^2 - 4 * a * c) - 1) / (2*a)
      dbeta3 = -(sqrt(b^2-4*a*c) - b)/(2*a^2) - c / (a * sqrt(b^2 - 4*a*c))
    }
    else if (solution == 2) {
      dbeta1 = 1/sqrt(b^2-4*a*c)
      dbeta2 = -(b/sqrt(b^2 - 4 * a * c) - 1) / (2*a)
      dbeta3 = c/(a*sqrt(b^2-4*a*c)) - (-sqrt(b^2 - 4*a*c))/(2*a^2)
    }
    else
      stop("Invalid solution")
    
    return(c(dbeta1, dbeta2, dbeta3))
  }
  else if (identical(powers, c(1, 2, 3))) { # cubic
    stop("Polynomial not yet supported")
  }
  else {
    stop("Polynomial not supported")
  }
  
}