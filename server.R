library(shiny)
library(waiter)

frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)

# load other source files
# and libraries
library(metaheuristicOpt)
library(DEoptim)
files.sources = list.files(path = "./R")
#cat(file=stderr(), files.sources)
files.sources = files.sources[files.sources != "app.R"] # don't call app
files.sources = paste("R/", files.sources, sep = "") # add directory name
sapply(files.sources, source)

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
  # x lower bound tries to be reasonably close to 0
  output$plot1 = renderPlot({
    ggp = ggplot(values$DT, aes(x = x, y = y)) +
      # geom_point(aes(color = color,
      #                shape = shape), size = 5) +
      geom_point(color = "red", shape = "circle", size = 5, alpha = 1) +
      lims(x = c(0.1, input$fp_bound), y = c(0, 1)) +
      theme_bw() + 
      # include so that colors don't change as more color/shape chosen
      #scale_color_discrete(drop = FALSE) +
      #scale_shape_discrete(drop = FALSE) +
      labs(y = "probability of response", x = "X",
           title = "Best FP fit")
    
    
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
    
    # make sure x and y are positive
    # this fixes missing values and is required for FP anyways
    x_coord = max(0.1, input$plot_click$x)
    y_coord = max(0, input$plot_click$y)
    add_row <- data.frame(x = x_coord,
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
    
    # fit either 2 or 3 degree polynomial
    if (input$fpdegree == 2) {
      out = fitted_logistic_fp2(successes, model_data$x, frac.powers)
    }
    else if (input$fpdegree == 3) {
      out = fitted_logistic_fp3(successes, model_data$x, frac.powers)
    }
    
    
    # save to reactive object
    values$DT$yhat = out$yhat
    
    # save data
    values$p1 = out$p1
    values$p2 = out$p2
    values$beta0 = out$beta0
    values$beta1 = out$beta1
    values$beta2 = out$beta2
    values$bound = input$fp_bound
    
    # save degree values
    if (input$fpdegree == 3) {
      values$beta3 = out$beta3
      values$p3  = out$p3
    }
    
  })
  
  # copy fitted model info to design inputs
  # needs something to show success
  observeEvent(input$copymodel, {
    updateNumericInput(session, "p1", value = values$p1)
    updateNumericInput(session, "p2", value = values$p2)
    updateNumericInput(session, "b0", value = values$beta0)
    updateNumericInput(session, "b1", value = values$beta1)
    updateNumericInput(session, "b2", value = values$beta2)
    updateNumericInput(session, "bound", value = values$bound)
    
    # set cubic options to NA if quadratic model is fit
    if (input$fpdegree == 2) {
      updateNumericInput(session, "p3", value = NA)
      updateNumericInput(session, "b3", value = NA)
    }
    else if (input$fpdegree == 3) {
      updateNumericInput(session, "p3", value = values$p3)
      updateNumericInput(session, "b3", value = values$beta3)
    }
    
    
  })
  
  output$model_out = renderPrint({
    
    # check for no model run
    if (length(values$beta0)==0) {
      # print("No model")
      cat("No Model\n")
    }
    else if (input$fpdegree == 2) {
      # print("p1:")
      # print(values$p1)
      # print("p2:")
      # print(values$p2)
      # print("beta0:")
      # print(values$beta0)
      # print("beta1:")
      # print(values$beta1)
      # print("beta2:")
      # print(values$beta2)
      cat("p1: ", values$p1, "\n",
          "p2: ", values$p2, "\n",
          "beta0: ", values$beta0, "\n",
          "beta1: ", values$beta1, "\n",
          "beta2: ", values$beta2, "\n",
          sep = ""
      )
    }
    else if (input$fpdegree == 3) {
      cat("p1: ", values$p1, "\n",
          "p2: ", values$p2, "\n",
          "p3: ", values$p3, "\n",
          "beta0: ", values$beta0, "\n",
          "beta1: ", values$beta1, "\n",
          "beta2: ", values$beta2, "\n",
          "beta3: ", values$beta3, "\n",
          sep = ""
      )
    }
  })
  
  
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
    # switch depending on if p3 or beta3 are missing
    if (is.na(input$p3) | is.na(input$b3)) {
      powers = as.numeric(c(input$p1, input$p2))
      betas = c(input$b0, input$b1, input$b2)
      degree = 2
    }
    else {
      powers = as.numeric(c(input$p1, input$p2, input$p3))
      betas = c(input$b0, input$b1, input$b2, input$b3)
      degree = 3
    }
    
    
    # algorithm options
    alg = metaheur_dict(input$alg)
    
    iter = input$iter
    swarm = input$swarm
    
    # design options
    pts = input$pts
    bound = input$bound
    
    # find optimal design
    od = ODpoly(powers, betas, alg, iter, swarm, pts, bound, degree)
    
    # store in reactive data
    values$OD$design = od$design
    values$OD$plot = od$plot
    values$OD$val = od$value
    
    
  })
  
  # update design output
  output$design_out = renderPrint({
    
    
    
    obj_val = values$OD$val
    raw = values$OD$design
    
    # case if algorithm hasn't run
    if (length(raw) == 0) {
      cat("No design") 
    }
    else { # all other cases
      
      # display objective value
      cat("-log(Det(M)) = ", obj_val, "\n", sep = "")
      #print(obj_val)
      
      l = length(raw)
      
      # purge points with zero weight
      if (sum(raw[(l/2 + 1):l]==0) > 0) {
        x_indices = which(raw[(l/2 + 1):l]==0)
        w_indices = x_indices + l/2
        raw = raw[,-c(x_indices, w_indices)]
        l = length(raw)
        cat("Purged points with weight 0\n")
      }
      
      
      # combine weights of identical points
      # sort as well
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
        #print("Combined identical points")
        cat("Combined identical points\n")
      }
      
      
      
      # labeling
      # probably a better way to do this
      # labs is a function => call it labbs
      labbs = character(l)
      for (i in 1:(l/2)) {
        labbs[i] = paste("x", toString(i), sep = "")
      }
      for (i in (l/2 + 1):l) {
        labbs[i] = paste("w", toString(i-l/2), sep = "")
      }
      
      # magic
      raw = c(raw)
      
      # sort by x's
      raw_x = raw[1:(l/2)]
      raw_w = raw[(l/2 + 1):l]
      r = rank(raw_x)
      raw_x = raw_x[r]
      raw_w = raw_w[r]
      raw = c(raw_x, raw_w)
      
      names(raw) = labbs
      
      out = raw
      #raw
      cat(labbs[1:(l/2)], "\n", sep = "    ")
      cat(round(out[1:(l/2)], 3), "\n")
      cat(labbs[(l/2 + 1):l], "\n", sep = "    ")
      cat(round(out[(l/2 + 1):l], 3))
    }
    
  })
}