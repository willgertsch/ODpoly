# testing for fractional polynomials
# get an idea of what they look like
powers = c(1/2,1/2)
betas = c(1, 1, 1)

xs = seq(0, 10, 0.1)

# compute x terms
if (powers[1] == powers[2]) { # apply correct x function if repeated power
  if (sum(xs<=0)>1) {
    print("Fractional polynomials are only defined on (0, inf]")
  }
  
  x1 = xs^powers[1]
  x2 = log(xs) * xs^powers[1]
} else {
  x1 = xs^powers[1]
  x2 = xs^powers[2]
}


# compute y
ys = betas[1] + betas[2] * x1 + betas[3] * x2

# plot
plot(ys~xs, type = "l")

for (i in seq(-1,1,0.1)) {
  betas = c(i, i, i)
  # compute x terms
  if (powers[1] == powers[2]) { # apply correct x function if repeated power
    if (sum(xs<=0)>1) {
      print("Fractional polynomials are only defined on (0, inf]")
    }
    
    x1 = xs^powers[1]
    x2 = log(xs) * xs^powers[1]
  } else {
    x1 = xs^powers[1]
    x2 = xs^powers[2]
  }
  
  
  # compute y
  ys = betas[1] + betas[2] * x1 + betas[3] * x2
  
  lines(ys ~ xs) # additional lines
}