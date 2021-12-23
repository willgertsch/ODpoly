# example usage of ODpoly
library(metaheuristicOpt)

set.seed(122021)

# model pararms
powers = c(2, -2)
betas = c(2, 1, -4)

# algorithm options
alg = "DE"
iter = 1000
swarm = 100

# design options
pts = 4
bound = 10

# find optimal design
od = ODpoly(powers, betas, alg, iter, swarm, pts, bound)

# design
od$design

# check equivalence theorem
od$plot
