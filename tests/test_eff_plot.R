# test_eff_plot.R
# test eff_plot function

# model pararms
powers = c(2, -2)
beta = c(2, 1, -4)

# design options
pts = 3
bound = 10

lam.grid = seq(0.0, 0.9, 0.1)

out = eff_plot(beta, powers, bound, pts, lam.grid, p = 0.5)
out$plot1
out$plot2
