library(shiny)
library(waiter)

# globals
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
metaheuristics = c("Particle Swarm Optimization", 
                   "Grey Wolf Optimizer", 
                   "Harmony Search Algorithm", 
                   "Moth Flame Optimizer",
                   "Differential Evolution")


ODpolyApp <- function(...) {
  ui <- fixedPage(
    withMathJax(),
    
    tabsetPanel(type = "pills", id = "tabpanel",
                
                tabPanel("Background", 
                         
                         tags$h1(
                           "Optimal Designs for Binary Dose-Response Experiments with Fractional Polynomials and Metaheuristics", 
                           style="text-align:center;"
                         ),
                         
                         tags$h3(
                           "Introduction",
                           style="text-align:center;"
                         ),
                         
                         tags$p("
           Binary dose-reponse experiments are common across many different fields such as toxicology, epidemiology, drug development, and enviromental science.
           The goal of these studies is to identify how the dosage relates to the probability of toxicity, efficacy, or some other event.
           For example, Price et. al. (1985) studied the effect of ethylene glycol on birth defects in mice. 
           "),
           
           tags$p("
           Because the shape of the dose-response curve is usually non-linear, polynomial models are often used to model the probability of the event.
           However, standard polynomials may not be flexible enough to model the true shape of the curve.
           Fractional polynomials are an extension of the standard polynomial that permits more flexibility by allowing the powers to vary.
           Therefore a logistic regression with a fractional polynomial predictor is a sensible choice for modeling the dose-response curve
           "),
           
           tags$p("
           Before running an experiment, it is important to think carefully about the design.
           A good design can avoid waste, high costs, and improve statistical efficiency.
           Given that we plan to analyze the data using a fractional polynomial logistic model, we can optimize the dose levels and the number of subjects at each dose level to plan the best experiment for that model.
           "),
           
           
           
           tags$p("
           This web app allows the user to find locally D-optimal designs for a degree 2 fractional polynomial logistic model using different metaheuristic optimization algorithms.
           The designs are locally optimal, so the user must suppply regression coefficients and powers based on prior knowledge of the true dose-response curve.
           These values may be obtained from literature or by inputing a curve shape in the app itself.
           "),
           
           tags$h3(
             "Fractional Polynomial Logistic Model",
             style="text-align:center;"
           ),
           tags$p("
    Royston and Altman (1994) proposed fractional polynomials as a more flexible generalization of standard polynomials.
    Let \\(X\\) be a positive predictor variable. 
    A fractional polynomial of degree \\( m \\) is defined as
    $$
    \\phi_m (X , \\mathbf{\\beta}, \\mathbf{p}) = \\sum^{m+1}_{j=1} \\beta_j H_j(X)
    $$
    where \\(\\mathbf{\\beta}\\) and \\(\\mathbf{p}\\) are vectors of regression coefficients and powers respectively. 
    The function \\(H_j(X)\\) is defined recursively as
    $$
    H_j(X) = \\begin{cases} 
    X^{(p_j)} & p_j \\neq p_{j-1}\\\\
    H_{j-1}(X) \\ln X & p_j = p_{j-1}
    \\end{cases}
    $$
    where \\(H_1(X)=1\\) and \\(p_1 = 0\\).
    The power of \\(X\\) in parentheses is shorthand for the Box-Tidwell transformation, which is defined as
    $$
    X^{(p_j)} = \\begin{cases} 
    X^{p_j} & p_j \\neq 0\\\\
    \\ln X & p_j = 0
    \\end{cases}
    $$
    Royston and Altman argue that \\(m=2\\) with the set of powers \\(\\mathcal{P} = \\{ -2, -1, -0.5, 0, 0.5, 1, 2, 3\\}\\) is sufficient for most applications. 
           "),
    
    tags$p("
    Suppose we have a binary outcome modeled as \\(y_i \\sim \\text{Bernoulli}(p_i)\\) for observations \\(i = 1, \\dots, n\\) and \\(0 \\leq p_i < 1\\). 
    We can add a predictor \\(x_i\\) into the model using the logit link function
    $$\\log\\left( \\frac{p_i}{1-p_i}\\right) = \\eta_i$$
    where \\(\\eta_i\\) is a fractional polynomial in \\(x_i\\). 
    Rearranging, we obtain an expression for the mean of \\(y_i\\) given \\(\\eta_i\\)
    $$
    p_i = \\frac{e^{\\eta_i}}{1+e^{\\eta_i}}= \\frac{1}{1+e^{-\\eta_i}} 
    $$
    The information matrix for this model may be expressed as
    $$
    M(\\beta, \\mathbf{p}) = \\sum_{i=1}^k p_i (1-p_i) f(x_i) f(x_i)' = X'WX
    $$
    where \\(W\\) is a diagonal weight matrix with entries \\(p_i (1-p_i)\\). 
    For a 2nd degree fractional polynomial, the design matrix \\(X\\) has rows \\(f(x_i)' = (1,H_2(x_i),H_3(x_i))\\)."
    ),
    
    tags$h3(
      "Optimal Design",
      style="text-align:center;"
    ),
    tags$p("
           An experimental design \\(\\xi\\) may be expressed as a collection of design points \\(x_1, \\dots, x_k\\) and weights \\(w_1, \\dots, w_k\\)  for a fixed sample size \\(N\\).
           Multiplying \\(w_i\\) by \\(N\\) gives the approximate number of samples at predictor value \\(x_i\\).
           The design \\(\\xi\\) is optimal if it maximizes or minimizes some function of the model information matrix \\(M(\\beta, \\mathbf{p})\\).
           For example, D-optimality minimizes
           $$
           \\Psi(M(\\beta, \\mathbf{p})) = -\\log (|M(\\beta, \\mathbf{p})|)
           $$
           The D-optimal design is the design that minimizes the volume of confidence ellipsoid for the regression parameters.
           Since the information matrix depends on values of \\(\\beta\\) and \\(\\mathbf{p}\\), the design is locally optimal.
           Values for \\(\\beta\\) and \\(\\mathbf{p}\\) can be chosen based on previous studies or theory.
           "),
    tags$p("
           To check if a design is locally D-optimal, we can use a result called the equivalence theorem. The theorem says that the design is optimal if
           $$
           ch(x) = \\frac{\\exp(\\eta)}{(1+\\exp(\\eta))^2} f(x)'M(\\beta, \\mathbf{p}) f(x) - p \\leq 0
           $$
           for all values of \\( x\\) in the design space with equality at the optimal design points and where \\(p\\) is the number of regression coefficients.
           Plotting \\(ch(x)\\) provides a graphical check of optimality.
           "),
    
    tags$h3(
      "Metaheurstic Optimization",
      style="text-align:center;"[]
    ),
    
    tags$p(
      "Recall that finding the optimal design requires maximizing or minimizing \\(\\Psi(M(\\beta, \\mathbf{p}))\\).
      Because the objective function for this model is complex, we use metaheuristic optimization algorithms.
      These algorithms are general purpose and are able to handle constraints, non-convexity, and a large number of variables.
      One downside is that there is no guarantee that these algorithms will converge to the optimal design, but the equivalence theorem provides an easy way to check if the optimum has been reached.
      "
    ),
    
    tags$p(
      "We have included several different metaheuristic algorithms for finding optimal designs.
      These include Particle Swarm Optimization, the Grey Wolf Optimizer, the Harmony Search Algorithm, the Moth Flame Optimizer, and Differential Evolution.
      The algorithms selected generally have good performance for this problem, but based on testing we recommend Differential Evolution"
    )
    
                ),
    tabPanel("Fractional Polynomials",
             
             titlePanel(
               "Fractional Polynomials on Probability Scale"
             ),
             
             sidebarLayout(
               sidebarPanel(
                 "Options",
                 numericInput("fp_bound", "Upper bound", 10, 1, NA, 1)
                 
               ),
               mainPanel(
                 plotOutput("plot1", click = "plot_click"),
                 actionButton("fit", "Fit"),
                 actionButton("rem_point", "Remove Last Point"),
                 actionButton("clear", "Clear all"),
                 verbatimTextOutput("model_out"),
                 actionButton("copymodel", "Copy model to design input")
               )
             )
             
             ),
    tabPanel("Design", id = "design",
             
             titlePanel(
               "Find the optimal design"
             ),
             sidebarLayout(
               sidebarPanel(
                 # inputs
                 selectInput("p1", "Power 1", frac.powers, selected = 2),
                 selectInput("p2", "Power 2", frac.powers, selected = -2),
                 selectInput("p3", "Power 3", c(frac.powers, NA), selected = NA),
                 numericInput("b0", "Beta0", 2, -Inf, Inf, 0.01), # bad idea to use inf? probably
                 numericInput("b1", "Beta1", 1, -Inf, Inf, 0.01),
                 numericInput("b2", "Beta2", -4, -Inf, Inf, 0.01),
                 numericInput("b3", "Beta3", NA, -Inf, Inf, 0.01),
                 selectInput("alg", "Algorithms", metaheuristics, selected = "Differential Evolution"),
                 numericInput("iter", "Iterations", 1000, 1, 10e7, 1),
                 numericInput("swarm", "Swarm size", 100, 1, 10e5, 1),
                 numericInput("pts", "Max design points", 4, 1, 10, 1),
                 numericInput("bound", "Upper bound", 10, 1, NA, 1)
               ),
               mainPanel(
                 #radioButtons("color", "Pick Color", c("Pink", "Green", "Blue")),
                 #selectInput("shape", "Select Shape:", c("Circle", "Triangle")),
                 #plotOutput("plot1", click = "plot_click"),
                 #actionButton("fit", "Fit"),
                 #actionButton("rem_point", "Remove Last Point"),
                 #actionButton("clear", "Clear all"),
                 #verbatimTextOutput("model_out"),
                 plotOutput("sens_plot"),
                 actionButton("find", "Find optimal design"),
                 waiter::use_waiter(),
                 verbatimTextOutput("design_out")
               )
             )
             )
    )
  )
  
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
        lims(x = c(1, input$fp_bound), y = c(0, 1)) +
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
  
      out = fitted_logistic_fp2(successes, model_data$x, frac.powers)
      
      # save to reactive object
      values$DT$yhat = out$yhat
      
      # save data
      values$p1 = out$p1
      values$p2 = out$p2
      values$beta0 = out$beta0
      values$beta1 = out$beta1
      values$beta2 = out$beta2
      values$bound = input$fp_bound
      
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
      
      # set cubic options to NA
      updateNumericInput(session, "p3", value = NA)
      updateNumericInput(session, "b3", value = NA)
      
    })
    
    output$model_out = renderPrint({
      
      # check for no model run
      if (length(values$beta0)==0) {
        print("No model")
      }
      else {
        print("p1:")
        print(values$p1)
        print("p2:")
        print(values$p2)
        print("beta0:")
        print(values$beta0)
        print("beta1:")
        print(values$beta1)
        print("beta2:")
        print(values$beta2)
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
      if (length(raw) == 0)
        out = "No design"
      else { # all other cases
        
        # display objective value
        print("-log(Det(M)) = ")
        print(obj_val)
        
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
        
        # sort by x's
        raw_x = raw[1:(l/2)]
        raw_w = raw[(l/2 + 1):l]
        r = rank(raw_x)
        raw_x = raw_x[r]
        raw_w = raw_w[r]
        raw = c(raw_x, raw_w)

        names(raw) = labs
        
        out = raw
      }
        
      out
    })
  }
  
  
  shinyApp(ui, server, ...)
}