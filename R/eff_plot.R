# eff_plot.R
# construct efficiency plots based on Cook and Wong (1993)
# returns efficiency plot and individual sensitivity plots and designs
# returns 2 efficiency plots
# plot1: plot the c objective value
# plot2: approximate efficiency by taking the best value of c objective as optimum
# beta: list of regression coef
# powers: list of fractional powers, skipping 0th power
# bound: upper bound for experimental data
# pts: number of design points
# lam.grid: grid of lambdas to use for plot
eff_plot = function(betas, powers, bound, pts,
                    lam.grid = seq(0.0, 1, 0.1), 
                    p = 0.5) {
  
  degree = length(betas) - 1
  
  # find optimal designs for each objective separately
  # D-optimal design
  alg = 'DE'
  iter = 1000
  swarm = 50
  crit = 'D'
  cat("Finding D optimal design...\n")
  out = ODpoly(powers, betas, alg, iter, swarm, pts, bound, degree, crit, p)
  Dopt_value = out$value
  Dplot = out$plot
  
  # EDp c-optimal design
  # optimal design is 1 design point at the value of EDp
  # not going to use efficiency because of matrix singularity
  # just use the objective value at each lambda
  
  # loop through specified lambda grid
  plots = list()
  designs = list()
  for (lam in lam.grid) {
    
    cat("Finding design for lambda = ", lam, "\n", sep = "")
    # find design for given lambda
    crit = "Dual"
    out = ODpoly(powers, betas, alg, iter, swarm, pts, bound, degree, crit, p,
                 lam = lam)
    # store result
    #append(plots, out$plot)
    #append(designs, out$design)
    plots[[length(plots) + 1]] = out$plot
    designs[[length(designs) + 1]] = out$design
  }
  
  
  # compute objective values for each design found
  # need to comppute objective function values
  D_func = obj_function_factory(powers, betas, degree, crit = 'D', bound, p, lam)
  C_func = obj_function_factory(powers, betas, degree, crit = 'EDp', bound, p, lam)
  
  # compute objective values for each design in list
  Dobj_values = unlist(lapply(designs, D_func))
  Cobj_values = unlist(lapply(designs, C_func))
  
  # data
  plot_data = data.frame(
    lambda = lam.grid,
    Dobj = Dobj_values,
    C_obj_value = exp(-Cobj_values)
  )
  plot_data$D_eff = (exp(plot_data$Dobj)/exp(Dopt_value))^(1/(degree + 1))
  
  
  # deal with ggplot's way of doing things
  min_c = min(plot_data$C_obj_value)
  plot_data2 = data.frame(
    lambda = c(lam.grid, lam.grid, lam.grid),
    value = c(plot_data$D_eff, plot_data$C_obj_value, min_c/plot_data$C_obj_value),
    type = c(rep("D-efficiency", length(lam.grid)),
             rep("var(c'beta)", length(lam.grid)),
             rep("Approx. c-efficiency", length(lam.grid)))
  )
  
  plt1 = ggplot(subset(plot_data2, type != "Approx. c-efficiency"), 
                aes(x = lambda, y= value, color = type)) +
    geom_point() + geom_line() +
    theme_bw() +
    labs(x = "lambda", y = "",
         title = "D-efficiency and C objective value by lambda") +
    theme(legend.title=element_blank()) +
    scale_x_continuous(breaks = lam.grid)
  
  plt2 = ggplot(subset(plot_data2, type != "var(c'beta)"), 
                aes(x = lambda, y= value, color = type)) +
    geom_point() + geom_line() +
    theme_bw() +
    labs(x = "lambda", y = "",
         title = "D-efficiency vs c-efficiency by lambda") +
    theme(legend.title=element_blank()) +
    scale_x_continuous(breaks = lam.grid)
  
  # return
  return(list(plot1 = plt1, plot2 = plt2, designs = designs, sens_plots = plots))
}