# plot sensitivity function
# xvals: vector of x values to plot
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
library(ggplot2)
plot_sens = function(xvals, vars, betas, powers) {
  
  # compute sens function
  yvals = sapply(xvals, sens, vars, betas, powers)
  
  # old base R plots
  # plot(yvals ~ xvals, type ="l")
  # abline(h=0)
  
  # replaced with ggplot because can return plot object
  p = ggplot(mapping = aes(y = yvals, x = xvals)) + 
    geom_line(color = "blue") + 
    geom_hline(yintercept = 0) +
    theme_bw() +
    labs(title = "Equivalence Theorem Check") +
    xlab("x") + ylab("ch(x)")
  
  return(p)
}