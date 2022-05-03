# eff_plot.R
# construct efficiency plots based on Cook and Wong (1993)
# beta: list of regression coef
# powers: list of fractional powers, skipping 0th power
# bound: upper bound for experimental data
# pts: number of design points
# lam.grid: grid of lambdas to use for plot
eff_plot = function(beta, powers, bound, pts,
                    lam.grid = seq(0.0, 1, 0.1), 
                    p = 0.5) {
  
  degree = length(beta) + 1
  
  # find optimal designs for each objective separately
  # D-optimal design
  alg = 'DE'
  iter = 1000
  swarm = 50
  crit = 'D'
  cat("Finding D optimal design...\n")
  out = ODpoly(powers, beta, alg, iter, swarm, pts, bound, degree, crit, p)
  Dcrit = out$value
  Dplot = out$plot
  
  # EDp c-optimal design
  # optimal design is 1 design point at the value of EDp
  # just use approximate solution with 2 points to avoid matrix singularity
  crit = 'EDp'
  out = ODpoly(powers, beta, alg, iter, swarm, pts = 2, bound, degree, crit, p)
  Ccrit = out$value
  Cplot = out$plot
  
  # loop through specified lambda grid
  obj_values = c()
  plots = list()
  designs = c()
  for (lam in lam.grid) {
    
    # find design for given lambda
    crit = "Dual"
    out = ODpoly(powers, beta, alg, iter, swarm, pts, bound, degree, crit, p,
                 lam = lam)
    # store result
    append(obj_values, out$value)
    append(plots, out$plot)
    append(designs, out$design)
  }
  
  # data
  result_data = data.frame(
    lambda = lam.grid,
    
  )
  
  # plot result
  main_plot = ggplot()
}