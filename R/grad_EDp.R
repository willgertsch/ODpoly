# grad_EDp.R
# finds EDp if it exists and then computes the gradient
# for a given fractional polynomial
# returns a list with EDp value and gradient
grad_EDp = function(beta, powers, bound, p) {
  
  force(powers)
  force(beta)
  
  
  # find EDp
  # define function to use for root finding
  # define a powers vector that includes 0
  zpowers = c(0, powers)
  if (length(beta) == 3) { # degree 2
    f = function(x) {
      beta1 = beta[1] - log(p/(1-p))
      x1 = H(2, x, zpowers)
      x2 = H(3, x, zpowers)
      beta1 + beta[2] * x1 + beta[3] * x2
    }
  }
  else if (length(beta) == 4) { # degree 3
    f = function(x) {
      beta1 = beta[1] - log(p/(1-p))
      x1 = H(2, x, zpowers)
      x2 = H(3, x, zpowers)
      x3 = H(4, x, zpowers)
      beta1 + beta[2] * x1 + beta[3] * x2 + beta[4] * x3
    }
  }
  else
    stop("Polynomial not supported for finding EDp")
  
  # find the roots
  library(rootSolve)
  roots = uniroot.all(f, c(0.1, bound), n=10000000, maxiter = 10000)
  
  # only interested in the minimum root found
  # return NA if no roots found
  if (length(roots) == 0)
    EDp = NA
  else
    EDp = min(roots)
  
  # compute gradient
  if (is.na(EDp))
    grad = NA
  else {
    
    if (length(beta) == 3) {
      Q = beta[2] * dH(2, EDp, zpowers) + beta[3] * dH(3, EDp, zpowers)
      
      # gradient calculation
      grad = numeric(3)
      grad[1] = -1/Q
      grad[2] = -H(2, EDp, zpowers)/Q
      grad[3] = -H(3, EDp, zpowers)/Q
    }
    else if (length(beta) == 4) {
      Q = beta[2] * dH(2, EDp, zpowers) + beta[3] * dH(3, EDp, zpowers) +
        beta[4] * dH(4, EDp, zpowers)
      
      grad = numeric(4)
      grad[1] = -1/Q
      grad[2] = -H(2, EDp, zpowers)/Q
      grad[3] = -H(3, EDp, zpowers)/Q
      grad[4] = -H(4, EDp, zpowers)/Q
    }
    else
      stop("Polynomial not supported for gradient computation")
  }
  
  # return
  return(list(EDp = EDp, grad = grad))
}