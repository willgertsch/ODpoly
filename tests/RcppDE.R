# RcppDE.R
# test of using this implementation of DE
# hopefully should be faster
library(RcppDE)

# degree 2
# model pararms
powers = c(2, -2)
betas = c(2, 1, -4)
degree = 2
obj_func = obj_function_factory(powers, betas, degree)
neg_obj_func = function(x) {
  return(-obj_func(x))
}

bound = 10
pts = 3

numVar = 2*pts
lower = c(rep(0.1, pts), rep(0, pts))
upper = c(rep(bound, pts), rep(1, pts))

control = RcppDE::DEoptim.control(NP = 100, itermax = 1000, trace = F)
set.seed(1234)
system.time({
  out = RcppDE::DEoptim(neg_obj_func, lower, upper, control)
})

out$optim

# cubic test
# parameters generated using the app
powers = c(1, 1, 1)
betas = c(3.192414, -7.733277, 7.02568, -1.62532)
degree = 3
obj_func = obj_function_factory(powers, betas, degree)

bound = 10
pts = 6

numVar = 2*pts
lower = c(rep(0.1, pts), rep(0, pts))
upper = c(rep(bound, pts), rep(1, pts))

control = RcppDE::DEoptim.control(NP = 200, itermax = 10000, trace = F)
set.seed(1234)
system.time({ # 77.8 seconds with 10000 iter
  out = RcppDE::DEoptim(neg_obj_func, lower, upper, control)
})

out$optim

# supposedly the Rccp package doesn't help much
# let's try the R only version
library(DEoptim)

# quadratic
powers = c(2, -2)
betas = c(2, 1, -4)
degree = 2
obj_func = obj_function_factory(powers, betas, degree)

bound = 10
pts = 3

numVar = 2*pts
lower = c(rep(0.1, pts), rep(0, pts))
upper = c(rep(bound, pts), rep(1, pts))

control = DEoptim::DEoptim.control(NP = 100, itermax = 1000, trace = F)
set.seed(1234)
system.time({
  out = DEoptim::DEoptim(neg_obj_func, lower, upper, control)
})

out$optim

# cubic
powers = c(1, 1, 1)
betas = c(3.192414, -7.733277, 7.02568, -1.62532)
degree = 3
obj_func = obj_function_factory(powers, betas, degree)

bound = 10
pts = 6

numVar = 2*pts
lower = c(rep(0.1, pts), rep(0, pts))
upper = c(rep(bound, pts), rep(1, pts))

control = DEoptim::DEoptim.control(NP = 200, itermax = 10000, trace = F)
set.seed(1234)
system.time({
  out = DEoptim::DEoptim(neg_obj_func, lower, upper, control)
}) # 84.4 seconds

out$optim

# conclusion: Rccp is a little faster, but maybe not worth it if init pop is borken
# trying an init population
# equally weighted random design
powers = c(1, 1, 1)
betas = c(3.192414, -7.733277, 7.02568, -1.62532)
degree = 3
obj_func = obj_function_factory(powers, betas, degree)

bound = 10
pts = 6

numVar = 2*pts
lower = c(rep(0.1, pts), rep(0, pts))
upper = c(rep(bound, pts), rep(1, pts))
init = matrix(rep(c(runif(pts, 0.1, bound),
                    rep(1/pts, pts)), swarm), nrow=swarm, byrow=T)

control = DEoptim::DEoptim.control(NP = 200, itermax = 1000, trace = F,
                                   initialpop = init)
set.seed(1234)
system.time({
  out = DEoptim::DEoptim(neg_obj_func, lower, upper, control)
})

# initial pop works! Will use this one instead