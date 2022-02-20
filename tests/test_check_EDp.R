# test_check_EDp.R
# testing out finding different EDps

# ED50
# (x-1)^2 = 0
# solution should be 1
p = 0.5
beta = c(1,-2, 1) 
powers = c(1, 2)
bound = 3
check_EDp(p, beta, powers, bound)

# change bound to exclude root
# should be NA
p = 0.5
beta = c(1,-2, 1) 
powers = c(1, 2)
bound = 0.2
check_EDp(p, beta, powers, bound)

# imaginary roots
# x^2 + x + 1 = 0
# should be NA
p = 0.5
beta = c(1,1, 1) 
powers = c(1, 2)
bound = 3
check_EDp(p, beta, powers, bound)

# root 0 and root of 5
# x(x-5) = 0
# should find 5
p = 0.5
beta = c(0, -5, 1) 
powers = c(1, 2)
bound = 10
check_EDp(p, beta, powers, bound)

# root 0 and root of -5
# x(x+5) = 0
# should be NA
p = 0.5
beta = c(0, 5, 1) 
powers = c(1, 2)
bound = 10
check_EDp(p, beta, powers, bound)

# changing percentiles
p = 0.9
beta = c(1,-2, 1) 
powers = c(1, 2)
bound = 3
check_EDp(p, beta, powers, bound)

