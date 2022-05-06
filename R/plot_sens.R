# plot sensitivity function
# xvals: vector of x values to plot
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
library(ggplot2)
plot_sens = function(xvals, vars, betas, powers, degree = 2, crit = "D", bound, p, lam) {
  
  # if c-optimal design, compute gradient
  if (crit == "EDp" | crit == "Dual") {
    EDp_grad = grad_EDp(betas, powers, bound, p = p)
    if (is.na(EDp_grad$EDp))
      stop("No X value for ED50 found.")
    else {
      dg = EDp_grad$grad
    }
  }
  
  # compute sens function
  yvals = sapply(xvals, sens, vars, betas, powers, degree, crit, bound, dg, lam)
  
  # old base R plots
  # plot(yvals ~ xvals, type ="l")
  # abline(h=0)
  
  # compute ch(x) at design points
  l = length(vars)
  design_points = vars[1:(l/2)]
  # pts_ch = sapply(design_points, sens, vars, betas, powers, degree)
  # actually better to set points at what the y-val is supposed to be
  pts_ch = rep(0, l/2)
  
  # replaced with ggplot because can return plot object
  p = ggplot(mapping = aes(y = yvals, x = xvals)) + 
    geom_line(color = "blue") + 
    geom_hline(yintercept = 0) +
    geom_point(aes(x = design_points, y = pts_ch), col = "red", size = 3) + 
    geom_vline(xintercept = design_points, color = "red", linetype = "dashed") +
    theme_bw() +
    labs(title = "Equivalence Theorem Check") +
    xlab("x") + ylab("ch(x)")
  
  return(p)
}