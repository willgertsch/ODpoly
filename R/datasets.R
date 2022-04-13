# datasets.R
# test datasets

# microbial risk data
# Black et al 1988
# Campylobacter Jejuni in Healthy Volunteers
# outcome is proportion infected
black = data.frame(
  y = c(5/10, 6/10, 11/13, 8/11, 15/19, 5/5, 4/4),
  x = c(8e2, 8e3, 9e4, 8e5, 1e6, 1e8, 1e8)
)

# Hornick et al 1970
# Salmonella Typhi Quailes in Healthy Human Subjects
# outcome is proportion infected
hornick = data.frame(
  y = c(0/14, 32/116, 16/32, 8/9, 40/42),
  x = c(1e3, 1e5, 1e7, 1e8, 1e9)
)

# Devopmental tox data
# Faes et al 2003
# outcome is proportion of malformations in live births
# shift dose by 0.1 to avoid model issues
# Ethylene glycol
EG = data.frame(
  y = c(.04, .667, .818, .957),
  x = c(0.1, 750.1, 1500.1, 3000.1)
)

# Diethylene glycol dimethyl ether
# dose shifted by 0.01
DYME = data.frame(
  y = c(.048, .05, .25, .826, 1),
  x = c(0.1, 62.6, 125.1, 250.1, 500.1)
)

# EPA report
# Weber CI, et al. 1989. Short-term methods for estimating the chronic toxicity of effluents and receiving waters to freshwater organisms, 2nd ed. EPA/600/4-89/001A. U.S. Environmental Pro-tection Agency, Cincinnati, OH.
# Survival of Fathead minnows exposed to NaPCP
# shift by 0.1
weber = data.frame(
  y = c(1, .8, .9, .9, .7, .4),
  x = c(0.1, 32.1, 64.1, 128.1, 256.1, 512.1)
)
