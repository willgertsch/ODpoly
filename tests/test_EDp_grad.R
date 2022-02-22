# test_EDp_grad.R
# test out gradient computation

# random quad from fp fitting app
# two possible ED50s
p = 0.5
beta = c(-3.118491, 1.963149, -0.1814034)
powers = c(1, 2)
bound = 10
check_EDp(p, beta, powers, bound)
EDp_grad(beta, powers, solution = 1)

# from test data
beta = c(-4.20494, 1.841736, -0.1489127)
check_EDp(p, beta, powers, bound) # should be around 2.5
EDp_grad(beta, powers, solution = 1)

# from article data
beta = c(-2.492062, 0.004352718, -8.517725e-07)
bound = 3000
check_EDp(p, beta, powers, bound) # around 600
EDp_grad(beta, powers, solution = 1)
