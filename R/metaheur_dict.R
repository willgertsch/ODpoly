# metaheur_dict.R
# dictionary function for metaheuristics
metaheur_dict = function(fullname){
  if (fullname == "Particle Swarm Optimization")
    return("PSO")
  else if (fullname == "Grey Wolf Optimizer")
    return("GWO")
  else if (fullname == "Harmony Search Algorithm")
    return("HS")
  else if (fullname == "Moth Flame Optimizer")
    return("MFO")
  else if (fullname == "Differential Evolution")
    return("DE")
  else if (fullname == "Fast DE")
    return("FDE")
}