# benchmark.R
# run comparison of algorithms and parameters
library(bench)

# parameter space
algorithms = c("DE", "PSO", "GWO", "HS")
iterations = c(50, 100, 200, 500, 1000)
swarms = c(20, 50, 100)
pts = c(3, 4, 5)

# problems
# problem 1: standard quad 3 pt design inspired from Fornius
beta1 = c(2, 0, -0.1)
power1 = c(1, 2)
bound1 = 10

# problem 2: standard quad 4 pt design inspired from Fornius
beta2 = c(-2, 10, -4)
power2 = c(1, 2)
bound2 = 5

# problem 3: standard cubic with most iconic shape
# obtained using plotting app
beta3 = c(1.456149, -0.2314059, 0.006198903, -4.310079e-05)
power3 = c(1, 2, 3)
bound3 = 100

# problem 4: FP hornick data
beta4 = c(0.1308073, 0.0001145817, -361.184)
power4 = c(0.5, -0.5)
bound4 = 1e9

# problem 5: FP DYME data, also has log term
beta5 = c(-5.142485, -0.8556875, 0.7279428)
power5 = c(0, 0.5)
bound5 = 500.1

# problem 6: FP cubic sigmoidal shape
beta6 = c(14.69976, -1217.967, -5154.445, 7129.894)
power6 = c(-1, -2, -2)
bound6 = 100

# test run
# out = ODpoly(power6, beta6, alg = "DE", iter = 750, swarm = 50, pts = 4, 
#        bound = bound6, degree = 3, crit = "D")
# out$plot

# storage
N = length(algorithms) * length(iterations) * length(swarms) * length(pts) * 6 * 5
sim_data = data.frame(
  algorithm = rep(NA, N),
  iterations = rep(NA, N),
  swarm_size = rep(NA, N),
  design_pts = rep(NA, N),
  problem = rep(NA, N),
  obj_value = rep(NA, N),
  time = rep(NA, N)
)

# start_time = Sys.time()
# out = ODpoly(power6, beta6, alg = "DE", iter = 1000, swarm = 100, pts = 4,
#        bound = bound6, degree = 3, crit = "D")
# end_time = Sys.time()
# end_time - start_time


# main loop
start_time = Sys.time()
set.seed(445)
i = 1
for (n in 1:5) {
  for (algo in algorithms) {
    for (iter in iterations) {
      for (ss in swarms) {
        for (pt in pts) {
          # problem 1
          start_t = Sys.time()
          out = ODpoly(power1, beta1, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound1, degree = 2, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          # problem 2
          start_t = Sys.time()
          out = ODpoly(power2, beta2, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound2, degree = 2, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          # problem 3
          start_t = Sys.time()
          out = ODpoly(power3, beta3, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound3, degree = 3, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          # problem 4
          start_t = Sys.time()
          out = ODpoly(power4, beta4, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound4, degree = 2, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          # problem 5
          start_t = Sys.time()
          out = ODpoly(power5, beta5, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound5, degree = 2, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          # problem 6
          start_t = Sys.time()
          out = ODpoly(power6, beta6, alg = algo, iter = iter, swarm = ss,
                            pts = pt, bound = bound6, degree = 3, crit = "D")
          val = out$value
          end_t = Sys.time()
          t = end_t - start_t
          
          sim_data[i, ] = c(algo, iter, ss, pt, 1, val, t)
          i = i + 1
          
          cat(i, "/", N, "\n")
        }
      }
    }
  }
}


end_time = Sys.time()
end_time - start_time

# process data
sim_data$iterations = as.numeric(sim_data$iterations)
sim_data$swarm_size = as.numeric(sim_data$swarm_size)
sim_data$design_pts = as.numeric(sim_data$design_pts)
sim_data$obj_value = as.numeric(sim_data$obj_value)
sim_data$time = as.numeric(sim_data$time)

write.csv(sim_data, "benchmarking/sim1.csv")
