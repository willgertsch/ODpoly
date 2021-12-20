# plot sensitivity function
# xvals: vector of x values to plot
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
plot_sens = function(xvals, vars, betas, powers) {
  
  # compute sens function
  yvals = sapply(xs, sens, vars, betas, powers)
  
  plot(yvals ~ xvals, type ="l")
  abline(h=0)
}