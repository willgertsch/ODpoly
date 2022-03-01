# fracpoly.R
# utility functions for fractional polynomials

# Box-Tidwell transformation
bt = function(X, p) {
  if (p != 0)
    return(X^p)
  else if (p == 0)
    return(suppressWarnings(log(X)))
} 

# derivative of Box-Tidwell
dbt = function(X, p) {
  if (p != 0)
    return(p * X^(p-1))
  else if (p == 0)
    return(1/X)
} 

# H function
# j: index
H = function(j, X, powers) {
  if (j == 1) # base case
    return(1)
  if (powers[j] != powers[j-1])
    return(bt(X, powers[j]))
  else if (powers[j] == powers[j-1])
    return(suppressWarnings(log(X)) * H(j-1, X, powers))
}

# derivative of the H function
dH = function(j, X, powers) {
  if (j == 1) # base case
    return(0)
  if (powers[j] != powers[j-1])
    return(dbt(X, powers[j]))
  else if (powers[j] == powers[j-1])
    return(suppressWarnings(log(X)) * dH(j-1, X, powers) + 
             H(j-1, X, powers)/X)
}

# calculates the fractional polynomial for given X, coefficients, powers
# m: degree
fracpoly = function(X, betas, powers, m) {
  
  y = 0
  
  for (j in 1:(m+1)) {
    y = y + betas[j] * H(j, X, powers)
  }
  
  return(y)
  
}