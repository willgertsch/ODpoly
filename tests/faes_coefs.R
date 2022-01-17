# faes_coefs.R
# create data from Faes et. al. (2003) and fit FP

# table 1 data
EG = data.frame(
  dose = c(0.1, 750, 1500, 3000), # shift to avoid log problems
  malformations = c(4, 67, 82, 96)
)

# table 2 data
DYME = data.frame(
  dose = c(0.1, 62.5, 125, 250, 500), # shift to avoid log problems
  malformations = c(5, 5, 25, 83, 100)
)

# fit models
powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
EG.mod = fitted_logistic_fp(EG$malformations, EG$dose, powers)
DYME.mod = fitted_logistic_fp(DYME$malformations, DYME$dose, powers)


# plots
par(mfrow=c(2,1))
plot(malformations ~ dose, data = EG,  ylab = "% malformations")
lines(EG.mod$yhat * 100 ~ EG$dose, col = "red")
title("EG data", xlab = "dose")

plot(malformations ~ dose, data = DYME, ylab = "% malformations")
lines(DYME.mod$yhat * 100 ~ DYME$dose, col = "red")
title("DYME data", xlab = "dose")
par(mfrow=c(1,1))

# print model information
EG.mod
DYME.mod



