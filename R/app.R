library(shiny)
library(waiter)

# globals
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
metaheuristics = c(
  "Differential Evolution",
  "Particle Swarm Optimization",
  "Grey Wolf Optimizer",
  "Harmony Search Algorithm"
  )


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
                           style="text-align:left;"
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
           Before running an experiment, it is important to think carefully about its design.
           A good design can avoid waste, minimize high costs, and improve statistical efficiency.
           Given that we plan to analyze the data using a fractional polynomial logistic model, we can optimize the dose levels and the number of subjects at each dose level to plan the best experiment for that model.
           "),
           
           
           
           tags$p("
           This web app allows the user to find locally D-optimal designs for degree 2 and 3 fractional polynomial logistic models using metaheuristic optimization algorithms.
           The designs are locally optimal, so the user must supply regression coefficients and powers based on prior knowledge of the true dose-response curve.
           These values may be obtained from literature or by using the model fitting tool included in the app.
           "),
           
           tags$h3(
             "Fractional Polynomial Logistic Model",
             style="text-align:left;"
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
    Some examples of degree 2 fractional polynomials on this set are:
           "),
    tags$ul(
      tags$li("\\( \\beta_1 + \\beta_2 X^{-1} + \\beta_3 X^{-2}\\)"),
      tags$li("\\( \\beta_1 + \\beta_2 \\sqrt{X} + \\beta_3 \\sqrt{X} \\ln (X)\\)"),
      tags$li("\\( \\beta_1 + \\beta_2 \\ln (X) + \\beta_3 (\\ln (X))^2\\)"
      )
    ),
    
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
      style="text-align:left;"
    ),
    tags$p("
           An experimental design \\(\\xi\\) may be expressed as a collection of design points \\(x_1, \\dots, x_k\\) and weights \\(w_1, \\dots, w_k\\)  for a fixed sample size \\(N\\).
           Multiplying \\(w_i\\) by \\(N\\) gives the approximate number of samples at predictor value \\(x_i\\).
           The design \\(\\xi\\) is optimal if it maximizes or minimizes some function of the model information matrix \\(M(\\beta, \\mathbf{p})\\).
           In this app, we use a dual objective design criterion defined as
           $$
           \\Psi_{dual}(M) = \\frac{(1-\\lambda)}{p} \\log (|M|) - \\lambda \\log \\left(\\nabla g(\\beta)' M^{-1} \\nabla g(\\beta) \\right)
           $$
           The dual objective optimizes estimation of a percentile of the dose response curve and general paramter estimation by combining c and D optimality criteria.
           The weight parameter \\(\\lambda\\) denotes relative importance of each objective and a D optimal design may be obtained by setting \\(\\lambda = 0\\).
           Since the information matrix depends on values of \\(\\beta\\) and \\(\\mathbf{p}\\), the design is locally optimal.
           Values for \\(\\beta\\) and \\(\\mathbf{p}\\) can be chosen based on previous studies or on other prior information.
           "),
    tags$p("
           To check if a design is locally optimal, we can use a result called the equivalence theorem. The theorem says that the design is optimal if
           $$
           ch(x) = \\frac{(1-\\lambda)}{p} \\sigma(x) b(x)' M^{-1} b(x) + \\lambda \\frac{\\left(\\sqrt{\\sigma(x)} b(x)' M^{-1} \\nabla g(\\beta)\\right)^2 }{\\nabla g(\\beta)' M^{-1} \\nabla g(\\beta)} - 1 \\leq 0
           $$
           for all values of \\( x\\) in the design space with equality at the optimal design points and where \\(p\\) is the number of regression coefficients.
           Plotting \\(ch(x)\\) provides a simple graphical check of optimality.
           "),
    
    tags$h3(
      "Metaheurstic Optimization",
      style="text-align:left;"[]
    ),
    
    tags$p(
      "Recall that finding the optimal design requires maximizing or minimizing an objective function.
      Because the objective function for this model is complex with many special cases, we use metaheuristic optimization algorithms.
      These algorithms are very flexible and are able to handle constraints, non-convexity, and a large number of variables.
      One downside is that there is no guarantee that these algorithms will converge to the optimal design, but the equivalence theorem mentioned earlier provides an easy way to check if the optimum has been reached.
      "
    ),
    
    tags$p(
      "We have included several different metaheuristic algorithms for finding optimal designs.
      These include Particle Swarm Optimization, the Grey Wolf Optimizer, the Harmony Search Algorithm, the Moth Flame Optimizer, and Differential Evolution.
      The algorithms selected generally have good performance for this problem, but based on testing we recommend Differential Evolution"
    )
    
                ),
    tabPanel("Polynomial Fitting",
             
             titlePanel(
               "Fit fractional polynomials"
             ),
             tags$p("
                    This app allows to user to fit fractional polynomials by either uploading a file or by entering data by clicking on the plot.  Usage is as follows:
                    "),
             tags$ol(
               tags$li("Select desired options for upper bound and degree of polynomial."),
               tags$li("Click on plot to generate data for the probability of response at each X value. 
                       It is also possible to upload a .csv file. This file should have columns named y and x with probabilities and dose levels respectively."),
               tags$li("Click the \"Fit\" button to fit a fractional polynomial model to the data. 
                       All fractional polynomials with powers from the set {2, -1, -1/2, 0, 1/2, 1, 2, 3} will be fit to the data and the model with the lowest AIC will be returned.
                       If either of the standard options are selected, a standard quadratic or cubic polynomial will be fit."
                       ),
               tags$li("If desired, clicking \"Copy model to design input\" will transfer the resulting model to the design tab to find the optimal design.")
               ),
             
             sidebarLayout(
               sidebarPanel(
                 "Options",
                 numericInput("fp_bound", "Upper bound", 10, 1, NA, 1),
                 selectInput("fpdegree", "Degree", c(2, 3, "Standard quadratic", "Standard cubic"), selected = 2)
                 
               ),
               mainPanel(
                 plotOutput("plot1", click = "plot_click"),
                 actionButton("fit", "Fit"),
                 actionButton("rem_point", "Remove Last Point"),
                 actionButton("clear", "Clear all"),
                 fileInput("upload", "Import data from file (see 2.)", accept = ".csv"),
                 verbatimTextOutput("model_out"),
                 actionButton("copymodel", "Copy model to design input")
               )
             )
             
             ),
    tabPanel("Design", id = "design",
             
             titlePanel(
               "Find the optimal design"
             ),
             tags$p(
               "This app allows the user to find dual objective designs for fractional polynomials using a variety of metaheuristic optimization algorithms.
               Usage is as follows:"
             ),
             tags$ol(
               tags$li("Set power and beta values from prior knowledge or use the fractional polynomial tab to generate values associated with a desired shape.
                       If power 3 and beta 3 are set to NA, then a degree 2 polynomial is assumed.
                       Otherwise, a cubic will be used."),
               tags$li("Choose algorithm, number of iterations, and swarm size. These may need to be tweaked based on the specific design problem.
                       If the default algorithm fails to converge, try increasing the number of iterations and/or swarm size.
                       Trying another algorithm may also improve performance on specific problems"),
               tags$li("Choose the maximum number of design points. If you specify more points than are in the optimal design, the algorithm will try to merge repeated points or those with low weights. 
                       Adding a few more points than are necessary can also help convergence to the optimal design."),
               tags$li("Select the desired upper bound for the design. This should be set based on the context of the problem. 
                       For example, in clinical trials the upper bound could be the maximum safe dosage.
                       Large upper bounds may cause issues, so transforming X to another scale maybe helpful"),
               tags$li("Select a percentile for the c-optimality criterion."),
               tags$li("Choose a lambda value to set the relative importance of each objective. Lambda = 0 denotes a D-optimal design while labmda = 1 denotes a c-optimal design."),
               tags$li("Click the \"Find\" button to find the optimal design given the inputs.
                       The algorithm should take 10-20 seconds to find the design for default algorithm options.
                       Design points and weights are displayed rounded to 3 decimal places"),
               tags$li("A plot of the ch(x) function will also be displayed. If ch(x) = 1 for all x, this indicates a matrix singularity and ch(x) cannot be displayed.")
               
             ),
             
             sidebarLayout(
               sidebarPanel(
                 # inputs
                 selectInput("p1", "Power 1", frac.powers, selected = 2),
                 selectInput("p2", "Power 2", frac.powers, selected = -2),
                 selectInput("p3", "Power 3", c(frac.powers, NA), selected = NA),
                 numericInput("b0", "Beta 0", 2, -Inf, Inf, 0.01), # bad idea to use inf? probably
                 numericInput("b1", "Beta 1", 1, -Inf, Inf, 0.01),
                 numericInput("b2", "Beta 2", -4, -Inf, Inf, 0.01),
                 numericInput("b3", "Beta 3", NA, -Inf, Inf, 0.01),
                 selectInput("alg", "Algorithms", metaheuristics, selected = "Differential Evolution"),
                 numericInput("iter", "Iterations", 500, 1, 10e7, 1),
                 numericInput("swarm", "Swarm size", 20, 1, 10e5, 1),
                 numericInput("pts", "Max design points", 3, 1, 10, 1),
                 numericInput("bound", "Upper bound", 10, 1, NA, 1),
                 #selectInput("crit", "Design Criterion", c("D", "EDp", "Dual"), selected = "D"),
                 numericInput("p", "ED Percentile", 0.5, 0.01, .99, 0.01),
                 numericInput("lam", "Lambda", 0.5, 0, 1, 0.01)
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
                 actionButton("find", "Find"),
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
            #"BIC: ", values$bic, "\n",
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
            #"BIC: ", values$bic, "\n",
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
  
  
  shinyApp(ui, server, ...)
}