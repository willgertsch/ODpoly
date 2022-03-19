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
  values$upload_state = NULL # flag for managing file uploads
  values$DT <- data.frame(x = numeric(),
                          y = numeric(),
                          yhat = numeric()
                          #color = factor(),
                          #shape = factor()
  )
  
  ## 2. Create a plot
  # x lower bound tries to be reasonably close to 0
  output$plot1 = renderPlot({
    
    # update data frame if there is user imported ata
    file_input()
    # update upper bound based on data
    if (length(values$DT$x > 0)) { # only correct if there is data
      if (input$fp_bound < max(values$DT$x)) {
        updateNumericInput(session, "fp_bound", value = max(values$DT$x))
      }
    }
    
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
           title = "Best FP fit") +
      geom_abline(slope = 0, intercept = 0.5, linetype = 4) +
      annotate(geom = "text", x = input$fp_bound, y = 0.55, label = "ED50")
    
    
    # if there non NA values for the predicted values, plot these as well
    if (sum(!is.na(values$DT$yhat)) > 0) {
      ggp = ggp + geom_line(aes(x=x, y=yhat))
    }
    
    
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
    else if (input$fpdegree == "Standard quadratic") {
      out = standard_quad(successes, model_data$x)
    }
    else if (input$fpdegree == "Standard cubic") {
      out = standard_cubic(successes, model_data$x)
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
    values$aic = out$aic
    values$bic = out$bic
    
    # save degree values
    if (input$fpdegree == 3 | input$fpdegree == "Standard cubic") {
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
    if (input$fpdegree == 2 | input$fpdegree == "Standard quadratic") {
      updateNumericInput(session, "p3", value = NA)
      updateNumericInput(session, "b3", value = NA)
    }
    else if (input$fpdegree == 3 | input$fpdegree == "Standard cubic") {
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
    else if (input$fpdegree == 2 | input$fpdegree == "Standard quadratic") {
      
      cat("p1: ", values$p1, "\n",
          "p2: ", values$p2, "\n",
          "beta0: ", values$beta0, "\n",
          "beta1: ", values$beta1, "\n",
          "beta2: ", values$beta2, "\n",
          "AIC: ", values$aic, "\n",
          "BIC: ", values$bic, "\n",
          sep = ""
      )
    }
    else if (input$fpdegree == 3 | input$fpdegree == "Standard cubic") {
      cat("p1: ", values$p1, "\n",
          "p2: ", values$p2, "\n",
          "p3: ", values$p3, "\n",
          "beta0: ", values$beta0, "\n",
          "beta1: ", values$beta1, "\n",
          "beta2: ", values$beta2, "\n",
          "beta3: ", values$beta3, "\n",
          "AIC: ", values$aic, "\n",
          "BIC: ", values$bic, "\n",
          sep = ""
      )
    }
  })
  
  # thanks to https://stackoverflow.com/a/44206615
  # plotting checks current upload state and updates data if there is user
  # submitted data available
  observeEvent(input$upload, {
    values$upload_state <- 'uploaded'
  })
  
  file_input <- reactive({
    if (is.null(values$upload_state)) {
      return(NULL)
    } else if (values$upload_state == 'uploaded') {
      
      # use this to update plot data
      import_data = check_import_data(read.csv(input$upload$datapath))
      values$DT = data.frame(
        y = import_data$y,
        x = import_data$x,
        yhat = rep(NA, nrow(import_data))
      )
      return(NULL)
    } 
  })
  
  
  
  ############################################################################
  # code for finding optimal design
  ############################################################################
  
  # set up reactive data structure
  values$OD <- list(
    design = numeric(),
    plot = ggplot(),
    msg = character()
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
    #crit = input$crit
    crit = "Dual" # always dual since D is special case
    alpha = input$alpha
    p = input$p
    lam = input$lam
    
    # to avoid crashing the app, need to check if EDp exists
    if (crit == "EDp" | crit == "Dual") {
      EDp_grad = grad_EDp(betas, powers, bound, p = p)
      
      if (is.na(EDp_grad$EDp)) {
        # do nothing
        values$OD$msg = "No X value found for EDp."
        values$OD$plot = ggplot()
        return()
      }
      # else continue
    }
    # find optimal design
    od = ODpoly(powers, betas, alg, iter, swarm, pts, bound, degree, crit, p, lam)
    
    # store in reactive data
    values$OD$msg = ""
    values$OD$design = od$design
    values$OD$plot = od$plot
    values$OD$val = od$value
  })
  
  # update design output
  output$design_out = renderPrint({
    
    obj_val = values$OD$val
    raw = values$OD$design
    # crit = input$crit
    crit = "Dual"
    p = input$p
    
    # case if algorithm hasn't run
    if (length(raw) == 0) {
      cat("No design") 
    }
    else if (values$OD$msg != "") {
      cat("No X value for EDp.")
    }
    else { # all other cases
      
      # display objective value
      if (crit == "D") {
        cat("log(Det(M)) = ", obj_val, "\n", sep = "")
      }
      else if (crit == "EDp") {
        cat("[EDp']^t M^{-1} EDp' = ", -obj_val, "\n", sep = "")
        # display value for EDp
        if (is.na(input$p3) | is.na(input$b3)) {
          powers = as.numeric(c(input$p1, input$p2))
          beta = c(input$b0, input$b1, input$b2)
        }
        else {
          powers = as.numeric(c(input$p1, input$p2, input$p3))
          beta = c(input$b0, input$b1, input$b2, input$b3)
        }
        
        ED50 = grad_EDp(beta, powers, input$bound, p = p)$EDp
        cat("EDp = ", ED50, "\n", sep = "")
        
      } else if (crit == "Dual") {
        cat("lambda*Cobj + (1-lambda)*Dobj = ", -obj_val, "\n", sep = "")
        # display value for EDp
        if (is.na(input$p3) | is.na(input$b3)) {
          powers = as.numeric(c(input$p1, input$p2))
          beta = c(input$b0, input$b1, input$b2)
        }
        else {
          powers = as.numeric(c(input$p1, input$p2, input$p3))
          beta = c(input$b0, input$b1, input$b2, input$b3)
        }
        
        ED50 = grad_EDp(beta, powers, input$bound, p)$EDp
        cat("EDp = ", ED50, "\n", sep = "")
      }
      
      l = length(raw)
      
      labbs = names(raw)
      cat(labbs[1:(l/2)], "\n", sep = "    ")
      cat(round(raw[1:(l/2)], 3), "\n")
      cat(labbs[(l/2 + 1):l], "\n", sep = "    ")
      cat(round(raw[(l/2 + 1):l], 3))
    }
    
  })
}
