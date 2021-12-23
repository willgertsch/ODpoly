# plot sensitivity function
# xvals: vector of x values to plot
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
plot_sens = function(xvals, vars, betas, powers) {
  
  # compute sens function
  yvals = sapply(xs, sens, vars, betas, powers)
  
  # old base R plots
  # plot(yvals ~ xvals, type ="l")
  # abline(h=0)
  
  # replaced with ggplot because can return plot object
  library(ggplot2)
  p = ggplot(mapping = aes(y = yvals, x = xvals)) + 
    geom_line(color = "blue") + 
    geom_hline(yintercept = 0) +
    theme_bw()
  
  return(p)
}