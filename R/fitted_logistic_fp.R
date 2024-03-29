# fitted_logistic_fp.R
# fits a logistic model with fractional polynomial link
# degree 2
# returns fitted values, coefs, and powers
# successes: number of successes out of 100
# x: predictor pattern corresponding to number of successes
# powers: set of powers to select model forms from
fitted_logistic_fp2 = function(successes, x, powers) {
  
  # loop over all values of fractional powers
  # find lowest AIC
  num_models = length(powers)^2
  result = c(1,1,1)
  for (p1 in powers) {
    for (p2 in powers) {
      
      # define powers vector
      zpowers = c(0, p1, p2)
      
      x1 = H(2, x, zpowers)
      x2 = H(3, x, zpowers)
      
      # fit model 
      mod.p1.p2 = glm(cbind(successes, 100 - successes) ~ x1 + x2,
                      family = binomial)
      
      # record AIC
      result = rbind(result, c(AIC(mod.p1.p2), p1, p2))
      
    }
  }
  
  # remove first row
  result = result[-1, ]
  
  # find min aic
  min_aic_index = which.min(result[, 1])
  
  # get every needed to refit model
  p1 = result[min_aic_index, 2]
  p2 = result[min_aic_index, 3]
  zpowers = c(0, p1, p2)
  
  x1 = H(2, x, zpowers)
  x2 = H(3, x, zpowers)
  
  # fit model 
  mod = glm(cbind(successes, 100 - successes) ~ x1 + x2,
            family = binomial)
  
  # add predicted values to plot
  # these are the predicted probabilities
  yhat = predict(mod, type = "response")
  
  
  beta = coef(mod)
  
  # return list
  out = list(
    yhat = yhat, 
    p1 = unname(p1),
    p2 = unname(p2),
    beta0 = unname(beta[1]),
    beta1 = unname(beta[2]),
    beta2 = unname(beta[3]),
    aic = AIC(mod),
    bic = BIC(mod)
    )
  
  return(out)
}

# degree 3 version
fitted_logistic_fp3 = function(successes, x, powers) {
  
  # loop over all values of fractional powers
  # find lowest AIC
  num_models = length(powers)^3
  result = c(1,1,1,1)
  for (p1 in powers) {
    for (p2 in powers) {
      for (p3 in powers) {
        # define powers vector
        zpowers = c(0, p1, p2, p3)
        
        x1 = H(2, x, zpowers)
        x2 = H(3, x, zpowers)
        x3 = H(4, x, zpowers)
        
        # fit model 
        mod.p1.p2.p3 = glm(cbind(successes, 100 - successes) ~ x1 + x2 + x3,
                        family = binomial)
        
        # record AIC
        result = rbind(result, c(AIC(mod.p1.p2.p3), p1, p2, p3))
      }
    }
  }
  
  # remove first row
  result = result[-1, ]
  
  # find min aic
  min_aic_index = which.min(result[, 1])
  
  # get every needed to refit model
  p1 = result[min_aic_index, 2]
  p2 = result[min_aic_index, 3]
  p3 = result[min_aic_index, 4]
  zpowers = c(0, p1, p2, p3)
  
  x1 = H(2, x, zpowers)
  x2 = H(3, x, zpowers)
  x3 = H(4, x, zpowers)
  
  # fit model 
  mod = glm(cbind(successes, 100 - successes) ~ x1 + x2 + x3,
            family = binomial)
  
  # add predicted values to plot
  # these are the predicted probabilities
  yhat = predict(mod, type = "response")
  
  
  beta = coef(mod)
  
  # return list
  out = list(
    yhat = yhat, 
    p1 = unname(p1),
    p2 = unname(p2),
    p3 = unname(p3),
    beta0 = unname(beta[1]),
    beta1 = unname(beta[2]),
    beta2 = unname(beta[3]),
    beta3 = unname(beta[4]),
    aic = AIC(mod),
    bic = BIC(mod)
  )
  
  return(out)
}

standard_quad= function(successes, x) {
  
  x1 = x
  x2 = x^2
  mod = glm(cbind(successes, 100 - successes) ~ x1 + x2,
            family = binomial)
  # add predicted values to plot
  # these are the predicted probabilities
  yhat = predict(mod, type = "response")
  
  
  beta = coef(mod)
  
  # return list
  out = list(
    yhat = yhat, 
    p1 = 1,
    p2 = 2,
    beta0 = unname(beta[1]),
    beta1 = unname(beta[2]),
    beta2 = unname(beta[3]),
    aic = AIC(mod),
    bic = BIC(mod)
  )
  
  return(out)
}

standard_cubic= function(successes, x) {
  
  x1 = x
  x2 = x^2
  x3 = x^3
  mod = glm(cbind(successes, 100 - successes) ~ x1 + x2 + x3,
            family = binomial)
  # add predicted values to plot
  # these are the predicted probabilities
  yhat = predict(mod, type = "response")
  
  
  beta = coef(mod)
  
  # return list
  out = list(
    yhat = yhat, 
    p1 = 1,
    p2 = 2,
    p3 = 3,
    beta0 = unname(beta[1]),
    beta1 = unname(beta[2]),
    beta2 = unname(beta[3]),
    beta3 = unname(beta[4]),
    aic = AIC(mod),
    bic = BIC(mod)
  )
  
  return(out)
}