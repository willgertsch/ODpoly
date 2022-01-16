# fracpoly_tests.R
# test to make sure degree 2 fractional polynomials are working correctly
# will generate all 36 possible models and then see if match output from function

X = seq(0.1, 10, 0.01)

powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

yvals = c() # storage

# generate response values
# included both permutations (p1, p2) and (p2, p1)
for (p1 in powers) {
  for (p2 in powers) {
    if (p1 != p2) { 
      if (p1 == 0)
        x1 = log(X)
      else
        x1 = X^p1
      
      if (p2 == 0)
        x2 = log(X)
      else
        x2 = X^p2
      
      y = 1 + 1 * x1 + 1 * x2 # fractional polynomial
      yvals = rbind(yvals, c(p1, p2, y))
    }
    else {
      if (p1 == 0) {
        x1 = log(X)
        x2 = log(X)^2
      }
      else {
        x1 = (X^p1)
        x2 = (X^p1) * log(X)
      }
      y = 1 + 1 * x1 + 1 * x2
      yvals = rbind(yvals, c(p1,p2,y))
    }
  }
}

# generate same data using fracpoly()
fracpoly_test = c()
betas = c(1,1,1)
for (p1 in powers) {
  for (p2 in powers) {
    y = fracpoly(X, betas, c(0, p1, p2), 2)
    fracpoly_test = rbind(yvals, c(p1, p2, y))
  }
}

# for some reason this duplicates the last row
# delete
fracpoly_test = fracpoly_test[-65,]

# fracpoly_test[,1:2]
# yvals[, 1:2]

# should both be 0
sum(fracpoly_test[,1] != yvals[, 1])
sum(fracpoly_test[,2] != yvals[, 2])

# check y results
# should be 0
sum(fracpoly_test != yvals)

