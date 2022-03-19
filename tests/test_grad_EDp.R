# test_grad_EDp.R

# x^2 - 8 + 15
# should have ED50 of 3 and grad = (0.5, 1.5, 4.5)
beta = c(15, -8, 1)
powers = c(1, 2)
bound = 10
p = 0.5
grad_EDp(beta, powers, bound, p)

# cubic
# (x-3)(x-5)^2 = x^3 - 13x^2 + 55x -75
# EDp should be 3
beta = c(-75, 55, -13, 1)
powers = c(1, 2, 3)
grad_EDp(beta, powers, bound, p)

# cubic with triple real root
# (x-1)(x-2)(x-3) = x^3 - 6x^2 +11x -6
beta = c(-6, 11, -6, 1)
grad_EDp(beta, powers, bound, p)

# no solutions
beta = c(1, 1, 1)
powers = c(1,2)
grad_EDp(beta, powers, bound, p)

# fractional polynomial testing ################################
# one possible root
# EDp around 2.5
beta = c(-1.362606, 1.459538, 0.00294904)
powers = c(0, 3)
grad_EDp(beta, powers, bound, p)

# quadratic shape with 2 possible values
# ED50 around 1.5
powers = c(-1, 3)
beta = c(-5.760749,9.899159,0.01067219)
grad_EDp(beta, powers, bound, p)

# repeated root with 2 possible values
# ED50 around 1
powers = c(0.5, 0.5)
beta = c(6.868488,-7.188789,2.555781)
grad_EDp(beta, powers, bound, p)

# no solution within bound
# well I tried to at least, says its around 0.45
powers = c(-0.5, 0.5)
beta = c(-13.13786,7.365509, 3.260565)
grad_EDp(beta, powers, bound, p)

# cubic with 3 real solutions
# ED50 around 1
powers = c(1, 2, 0.5)
betas = c(-12.22262, -7.951942, 0.3482211, 19.17816)
grad_EDp(beta, powers, bound, p)


# logistic shape
# triple repeated power
# ED50 around 5.3 => says around .47
powers = c(-2, -2, -2)
beta = c(7.527883, 5.968856, -33.67614, -58.62805)
grad_EDp(beta, powers, bound, p)

# no solution within bound
powers = c(2, 2, 2)
beta = c(1.002221, 0.545415, -0.601232, 0.1635648)
grad_EDp(beta, powers, bound, p)

# problematic real world issue
powers = c(1,2,3)
beta = c(-23.1, 17.7, -4.7, 0.41)
bound = 5
p = 0.25
grad_EDp(beta, powers, bound, p)

