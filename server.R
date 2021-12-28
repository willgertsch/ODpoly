library(shiny)
library(waiter)

frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
xs = seq(0, 10, 0.01)


library(metaheuristicOpt)
ODpoly = function(powers, betas, alg = "DE", iter, swarm, pts, bound) {
  
  # define objective function
  obj_func = obj_function_factory(powers, betas)
  
  # set up for metaheuristics
  numVar = 2*pts # each point has one weight and last weight is sum of the others
  d = c(rep(c(1, bound), pts), rep(c(0,1), pts))
  rangeVar = matrix(d, nrow = 2)
  
  # algorithm params
  control = list(numPopulation = swarm, maxIter = iter)
  
  # find optimal design
  sol = metaOpt(obj_func, optimType = "MAX", algorithm = alg, 
                numVar, rangeVar, control)
  
  # save output
  result = sol$result
  
  # plot sensitivity function
  xs = seq(1, bound, 0.01)
  p = plot_sens(xs, sol$result, betas, powers)
  
  # return
  out = list(design = result, plot = p)
}

obj_function_factory = function(powers, betas) {
  
  # check input
  if (length(powers) != length(betas)-1) # make sure there is a coefficient for each power
    return(0)
  
  force(powers)
  force(betas)
  
  lbeta = length(betas)
  
  # construct objective function
  # only input should be design on that point
  # vars is a list with the current x values and weights
  obj_func = function(vars) {
    
    # distinguish between points and weights
    pts = length(vars)/2
    x = vars[1:pts]
    w = vars[(pts+1):(2*pts)]
    s = sum(w, na.rm = T) # needed to fix if statement error
    if (s < 0 | s > 1) # constraint implementation
      return(-Inf)
    
    
    # compute x terms
    if (powers[1] == powers[2]) { # apply correct x function if repeated power
      if (sum(x<=1)>0) {
        #print("Fractional polynomials are only defined on (0, inf]")
        x[x<=1] = 1 # put back in design interval
      }
      
      x1 = x^powers[1]
      x2 = log(x) * x^powers[1]
    }
    else {
      x1 = x^powers[1]
      x2 = x^powers[2]
    }
    
    
    # compute eta
    eta = betas[1] + betas[2] * x1 + betas[3] * x2
    
    
    # weight function
    sigma = exp(eta)/(1+exp(eta))^2
    
    # information matrix
    # currently quadratic
    M = 0
    for (i in 1:pts) {
      
      # will need to update to use correct x functions
      m12 = x1[i]
      m13 = x2[i]
      m23 = x1[i]*x2[i]
      
      M_i = w[i] * sigma[i] * matrix(c(
        1, m12, m13,
        m12, x1[i]^2, m23,
        m13, m23, x2[i]^2
      ), ncol=3)
      
      M = M + M_i
      
    }
    
    
    # use information matrix to compute objective
    # using D for now
    # silence warnings
    
    obj_value = suppressWarnings(log(det(M)))
    
    
    
    # deal with NAs
    if (is.na(obj_value)) {
      return(-Inf)
    }
    else
      return(obj_value)
  }
  
  return(obj_func)
}

sens = function(z, vars, betas, powers) {
  
  # distinguish between points and weights
  pts = length(vars)/2
  x = vars[1:pts]
  w = vars[(pts+1):(2*pts)]
  s = sum(w, na.rm = T) # need to get rid of NAs here
  if (s < 0 | s > 1) # constraint implementation
    return(-Inf)
  
  # compute x and z terms using correct functional
  if (powers[1] == powers[2]) { # apply correct x function if repeated power
    # not sure I need this
    if (sum(x<=1)>0) {
      x[x<=1] = 1 # put back in design interval
    }
    
    x1 = x^powers[1]
    x2 = log(x) * x^powers[1]
    z1 = z^powers[1]
    z2 = log(z) * z^powers[2]
  }
  else {
    x1 = x^powers[1]
    x2 = x^powers[2]
    z1 = z^powers[1]
    z2 = z^powers[2]
  }
  
  # compute eta
  # additional log terms not implemented yet
  eta = betas[1] + betas[2] * x1 + betas[3] * x2
  etaz = betas[1] + betas[2] * z1 + betas[3] * z2
  
  # weight functions
  sigma = exp(eta)/(1+exp(eta))^2
  sigmaz = exp(etaz)/(1+exp(etaz))^2
  
  # compute information matrix
  # information matrix
  # currently quadratic
  M = 0
  for (i in 1:pts) {
    
    # will need to update to use correct x functions
    m12 = x1[i]
    m13 = x2[i]
    m23 = x1[i] * x2[i]
    
    M_i = w[i] * sigma[i] * matrix(c(
      1, m12, m13,
      m12, x1[i]^2, m23,
      m13, m23, x2[i]^2
    ), ncol=3)
    
    M = M + M_i
    
  }
  
  # compute matrix inverse and then sensitivity function
  # avoid singular matrices
  # solution from https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
  if (class(try(solve(M),silent=T))[1]!="matrix") {
    # set Minv to something
    #Minv = matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)
    y = 1
  }
  else {
    # compute sensitivity function
    Minv = solve(M)
    b = c(1, z1, z2)
    y = sigmaz * t(b) %*% Minv %*% b - 3
  }
  
  
  
  return(y)
  
  
}

# plot sensitivity function
# xvals: vector of x values to plot
# vars: solution vector
# betas: vector of coefficients
# powers: vector of powers
library(ggplot2)
plot_sens = function(xvals, vars, betas, powers) {
  
  # compute sens function
  yvals = sapply(xvals, sens, vars, betas, powers)
  
  # old base R plots
  # plot(yvals ~ xvals, type ="l")
  # abline(h=0)
  
  # replaced with ggplot because can return plot object
  p = ggplot(mapping = aes(y = yvals, x = xvals)) + 
    geom_line(color = "blue") + 
    geom_hline(yintercept = 0) +
    theme_bw() +
    labs(title = "Equivalence Theorem Check") +
    xlab("x") + ylab("Standardized variance")
  
  return(p)
}

server <- function(input, output, session) {
  
  # code for fractional polynomial plot
  # thanks to this person: http://www.ostack.cn/?qa=733734/
  ## 1. set up reactive dataframe
  values <- reactiveValues()
  values$DT <- data.frame(x = numeric(),
                          y = numeric(),
                          yhat = numeric()
                          #color = factor(),
                          #shape = factor()
  )
  
  ## 2. Create a plot
  output$plot1 = renderPlot({
    ggp = ggplot(values$DT, aes(x = x, y = y)) +
      # geom_point(aes(color = color,
      #                shape = shape), size = 5) +
      geom_point(color = "red", shape = "circle", size = 5, alpha = 1) +
      lims(x = c(1, 10), y = c(0, 1)) +
      theme_bw() + 
      # include so that colors don't change as more color/shape chosen
      #scale_color_discrete(drop = FALSE) +
      #scale_shape_discrete(drop = FALSE) +
      labs(y = "probability of response", x = "dose",
           title = "Click plot to enter data.\nClick \"Fit\" to fit a fractional polynomial to the data")
    
    
    # if there non NA values for the predicted values, plot these as well
    if (sum(!is.na(values$DT$yhat)) > 0)
      ggp = ggp + geom_line(aes(x=x, y=yhat))
    
    # display plot
    ggp
  })
  
  ## 3. add new row to reactive dataframe upon clicking plot ##
  observeEvent(input$plot_click, {
    # each input is a factor so levels are consistent for plotting characteristics
    # I modify this code to produce an outcome for logistic regression
    # each point represents 100 observations at that x value
    # y is the number of positive responses received
    # this gives a binomial response
    
    # make sure y is positive
    y_coord = max(0, input$plot_click$y)
    add_row <- data.frame(x = input$plot_click$x,
                          y = y_coord,
                          yhat = NA
                          #color = factor(input$color, levels = c("Pink", "Green", "Blue")),
                          #shape = factor(input$shape, levels = c("Circle", "Triangle"))
    )
    # add row to the data.frame
    values$DT <- rbind(values$DT, add_row)
  })
  
  ## 4. remove row on actionButton click ##
  observeEvent(input$rem_point, {
    rem_row <- values$DT[-nrow(values$DT), ]
    values$DT <- rem_row
  })
  
  # clear data frame
  observeEvent(input$clear, {
    values$DT <- data.frame(x = numeric(),
                            y = numeric(),
                            yhat = numeric()
                            #color = factor(),
                            #shape = factor()
    )
  })
  
  # find fractional polynomial algorithm
  observeEvent(input$fit, {
    
    # save model data
    model_data = values$DT
    
    # calculate number of successes
    successes = round(model_data$y * 100)
    
    # loop over all values of fractional powers
    # find lowest AIC
    num_models = length(frac.powers)^2
    result = c(1,1,1)
    for (p1 in frac.powers) {
      for (p2 in frac.powers) {
        
        # create x variables
        x1 = model_data$x^p1
        if (p1 == p2)
          x2 = log(model_data$x) * model_data$x^p2
        else
          x2 = model_data$x^p2
        
        
        
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
    
    # create x variables
    x1 = model_data$x^p1
    if (p1 == p2)
      x2 = log(model_data$x) * model_data$x^p2
    else
      x2 = model_data$x^p2
    
    # fit model 
    mod = glm(cbind(successes, 100 - successes) ~ x1 + x2,
              family = binomial)
    
    # add predicted values to plot
    # these are the predicted probabilities
    yhat = predict(mod, type = "response")
    
    # save to reactive object
    values$DT$yhat = yhat
    
    
    # return
    # give powers labels
    # out = c(coef(mod), p1, p2)
    # names(out) = c(names(out)[1:length(coef(mod))], 'p1', 'p2')
    # out
    
    # change values in inputs
    # have to get rid of names
    p1 = unname(p1)
    p2 = unname(p2)
    beta = coef(mod)
    beta0 = unname(beta[1])
    beta1 = unname(beta[2])
    beta2 = unname(beta[3])
    updateNumericInput(session, "p1", value = p1)
    updateNumericInput(session, "p2", value = p2)
    updateNumericInput(session, "b0", value = beta0)
    updateNumericInput(session, "b1", value = beta1)
    updateNumericInput(session, "b2", value = beta2)
  })
  
  # output$model_out <- renderPrint(
  #   mod()
  # )
  
  ############################################################################
  # code for finding optimal design
  ############################################################################
  
  # set up reactive data structure
  values$OD <- list(
    design = numeric(),
    plot = ggplot()
  )
  
  output$sens_plot = renderPlot({
    
    # load from reactive data
    ggp = values$OD$plot
    
    # display plot
    ggp
  })
  
  # action for find button
  observeEvent(input$find, {
    
    # set up spinny thing
    waiter <- waiter::Waiter$new(id = "sens_plot",
                                 html = spin_terminal(),
                                 color = "grey"
    )$show()
    waiter$show()
    on.exit(waiter$hide())
    
    # model pararms
    powers = as.numeric(c(input$p1, input$p2))
    betas = c(input$b0, input$b1, input$b2)
    
    # algorithm options
    alg = input$alg
    iter = input$iter
    swarm = input$swarm
    
    # design options
    pts = input$pts
    bound = input$bound
    
    # find optimal design
    od = ODpoly(powers, betas, alg, iter, swarm, pts, bound)
    
    # store in reactive data
    values$OD$design = od$design
    values$OD$plot = od$plot
    
    
  })
  
  # update design output
  output$design_out = renderPrint({
    
    raw = values$OD$design
    
    # case if algorithm hasn't run
    if (length(raw) == 0)
      out = "No design"
    else { # all other cases
      
      l = length(raw)
      
      # purge points with zero weight
      if (sum(raw[(l/2 + 1):l]==0) > 0) {
        x_indices = which(raw[(l/2 + 1):l]==0)
        w_indices = x_indices + l/2
        raw = raw[,-c(x_indices, w_indices)]
        l = length(raw)
        print("Purged points with weight 0")
      }
      
      
      # combine weights of identical points
      xs = raw[1:(l/2)]
      ws = raw[(l/2+1):l]
      if (length(unique(xs)) != length(xs)) {
        dups = xs[duplicated(xs)] # keep track of dups
        for (d in dups) { # iterate through and combine weights
          indices = xs == d
          new_w = sum(ws[indices]) # add w's for a specific duplicate
          ws[indices] = new_w # update w's; will drop later
        }
        xs = unique(xs) # remove duplicates
        ws = unique(ws)
        
        raw = c(xs, ws)
        l = length(raw)
        print("Combined identical points")
      }
      
      
      
      # labeling
      # probably a better way to do this
      labs = character(l)
      for (i in 1:(l/2)) {
        labs[i] = paste("x", toString(i), sep = "")
      }
      for (i in (l/2 + 1):l) {
        labs[i] = paste("w", toString(i-l/2), sep = "")
      }
      
      # magic
      raw = c(raw)
      
      names(raw) = labs
      
      
      out = raw
    }
    
    out
  })
}