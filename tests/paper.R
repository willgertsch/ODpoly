# paper.R
# run all the analyses seen in the paper

################################################################################
# things I need for tables
# 3 designs: D-optimal, dual, original
# d and approx c-efficiencies for both
# will start with the original design and return table of optimal design
# input raw data data frame with y = prob response and x = dose, n = number

# weber data ###################################################################
original_data = data.frame(
  y = weber$y,
  x = weber$x,
  n = rep(40, 6)
)

frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)

# convert to weighted design format
N = sum(original_data$n)
original_data$w = original_data$n/N
design = c(original_data$x, original_data$w)

# fit model to data and get parameter estimates
mod = fitted_logistic_fp2(round(100 * original_data$y), original_data$x, frac.powers)
powers = c(mod$p1, mod$p2)
betas = c(mod$beta0, mod$beta1, mod$beta2)
iter = 1000
swarm = 50
pts = 3
bound = max(original_data$x)

# evaluate for dual objective on a lambda grid
# also computes the D optimal design as well
eff_out = eff_plot(betas, powers, bound, pts,
                   lam.grid = seq(0.0, .9, 0.1), 
                   p = 0.5,
                   alg = "DE")

# compute design efficiency for original design
# create objective functions
D_func = obj_function_factory(powers, betas, crit = 'D', bound)
C_func = obj_function_factory(powers, betas, crit = 'EDp', bound, p=0.5)

# evaluate objective value
Dobj_val = D_func(design)
Cobj_val = exp(-C_func(design))

# best values obtained in grid evaluation
Doptimal = max(eff_out$obj_vals$Dobj)
Coptimal = min(eff_out$obj_vals$C_obj_value)


# compute efficiencies
Deff = (exp(Dobj_val)/exp(Doptimal))^(1/3) # for a degree 2 model
Ceff = Coptimal/Cobj_val

# return
original_design = c(design, Deff, Ceff)
# D optimal design is first stored in array
Doptimal_design = c(eff_out$designs[[1]], 1, Coptimal/eff_out$obj_vals$C_obj_value[1])
# get dual designs by looking at eff_out$dat
eff_out$dat
# lambda = 0.5 has Deff = 85% and Ceff = 75%
eff_out$designs[[6]]


# DYME #########################################################################
original_data = data.frame(
  y = DYME$y,
  x = DYME$x,
  n = c(21, 20, 24, 23, 23)
)

frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)

# convert to weighted design format
N = sum(original_data$n)
original_data$w = original_data$n/N
design = c(original_data$x, original_data$w)

# fit model to data and get parameter estimates
mod = fitted_logistic_fp2(round(100 * original_data$y), original_data$x, frac.powers)
powers = c(mod$p1, mod$p2)
betas = c(mod$beta0, mod$beta1, mod$beta2)
iter = 1000
swarm = 50
pts = 3
bound = max(original_data$x)

# evaluate for dual objective on a lambda grid
# also computes the D optimal design as well
eff_out = eff_plot(betas, powers, bound, pts,
                   lam.grid = seq(0.0, .9, 0.1), 
                   p = 0.5,
                   alg = "DE")

# compute design efficiency for original design
# create objective functions
D_func = obj_function_factory(powers, betas, crit = 'D', bound)
C_func = obj_function_factory(powers, betas, crit = 'EDp', bound, p=0.5)

# evaluate objective value
Dobj_val = D_func(design)
Cobj_val = exp(-C_func(design))

# best values obtained in grid evaluation
Doptimal = max(eff_out$obj_vals$Dobj)
Coptimal = min(eff_out$obj_vals$C_obj_value)

# compute efficiencies
Deff = (exp(Dobj_val)/exp(Doptimal))^(1/3) # for a degree 2 model
Ceff = Coptimal/Cobj_val

# return
original_design = c(design, Deff, Ceff)
# D optimal design is first stored in array
Doptimal_design = c(eff_out$designs[[1]], 1, Coptimal/eff_out$obj_vals$C_obj_value[1])
# get dual designs by looking at eff_out$dat
eff_out$dat
# lambda = 0.5 has Deff = 82% and Ceff = 73%
eff_out$designs[[6]]

# Hornick ######################################################################
original_data = data.frame(
  y = hornick$y,
  x = log10(hornick$x),
  n = c(14, 116, 32, 9, 42)
)

frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)

# convert to weighted design format
N = sum(original_data$n)
original_data$w = original_data$n/N
design = c(original_data$x, original_data$w)

# fit model to data and get parameter estimates
mod = fitted_logistic_fp2(round(100 * original_data$y), original_data$x, frac.powers)
powers = c(mod$p1, mod$p2)
betas = c(mod$beta0, mod$beta1, mod$beta2)
iter = 1000
swarm = 50
pts = 3
bound = max(original_data$x)

# evaluate for dual objective on a lambda grid
# also computes the D optimal design as well
eff_out = eff_plot(betas, powers, bound, pts,
                   lam.grid = seq(0.0, .9, 0.1), 
                   p = 0.5,
                   alg = "DE")

# compute design efficiency for original design
# create objective functions
D_func = obj_function_factory(powers, betas, crit = 'D', bound)
C_func = obj_function_factory(powers, betas, crit = 'EDp', bound, p=0.5)

# evaluate objective value
Dobj_val = D_func(design)
Cobj_val = exp(-C_func(design))

# best values obtained in grid evaluation
Doptimal = max(eff_out$obj_vals$Dobj)
Coptimal = min(eff_out$obj_vals$C_obj_value)

# compute efficiencies
Deff = (exp(Dobj_val)/exp(Doptimal))^(1/3) # for a degree 2 model
Ceff = Coptimal/Cobj_val

# return
original_design = c(design, Deff, Ceff)
# D optimal design is first stored in array
Doptimal_design = c(eff_out$designs[[1]], 1, Coptimal/eff_out$obj_vals$C_obj_value[1])
# get dual designs by looking at eff_out$dat
eff_out$dat
# lambda = 0.5 has Deff = 80% and Ceff = 72%
eff_out$designs[[6]]
