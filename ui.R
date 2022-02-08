library(shiny)
library(waiter)

# globals
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
metaheuristics = c(
  "Differential Evolution",
  "Particle Swarm Optimization",
  "Fast DE",
  "Grey Wolf Optimizer",
  "Harmony Search Algorithm",
  "Moth Flame Optimizer"
)

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
           Before running an experiment, it is important to think carefully about its design.
           A good design can avoid waste, minimize high costs, and improve statistical efficiency.
           Given that we plan to analyze the data using a fractional polynomial logistic model, we can optimize the dose levels and the number of subjects at each dose level to plan the best experiment for that model.
           "),
           
           
           
           tags$p("
           This web app allows the user to find locally D-optimal designs for degree 2 and 3 fractional polynomial logistic models using metaheuristic optimization algorithms.
           The designs are locally optimal, so the user must suppply regression coefficients and powers based on prior knowledge of the true dose-response curve.
           These values may be obtained from literature or by using the model fitting tool included in the app.
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
           Values for \\(\\beta\\) and \\(\\mathbf{p}\\) can be chosen based on previous studies or on other prior information.
           "),
    tags$p("
           To check if a design is locally D-optimal, we can use a result called the equivalence theorem. The theorem says that the design is optimal if
           $$
           ch(x) = \\frac{\\exp(\\eta)}{(1+\\exp(\\eta))^2} f(x)'M(\\beta, \\mathbf{p}) f(x) - p \\leq 0
           $$
           for all values of \\( x\\) in the design space with equality at the optimal design points and where \\(p\\) is the number of regression coefficients.
           Plotting \\(ch(x)\\) provides a simple graphical check of optimality.
           "),
    
    tags$h3(
      "Metaheurstic Optimization",
      style="text-align:center;"[]
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
    tabPanel("Fractional Polynomials",
             
             titlePanel(
               "Fractional Polynomials on Probability Scale"
             ),
             
             tags$p("
                    The coefficients and powers for the fractional polynomial model can be difficult to interpret, especially in the case of modeling the response probablity.
                    The probability of response for the logistic model with fractional polynomial predictor \\( \\eta \\) is
                    $$
                    p = \\frac{1}{1 + e^{-\\eta}}
                    $$
                    Thus it is difficult to match parameter values to the desired shape of the response curve.
                    "),
             tags$p("
                    This app allows the user to obtain fractional polynomial parameters for a desired curve shape. Usage is as follows:
                    "),
             tags$ol(
               tags$li("Select desired options for upper bound and degree of polynomial."),
               tags$li("Click on plot to generate data for the probability of response at each X value."),
               tags$li("Click the \"Fit\" button to fit a fractional polynomial model to the data. 
                       All fractional polynomials with powers from the set {2, -1, -1/2, 0, 1/2, 1, 2, 3} will be fit to the data and the model with the lowest AIC will be returned"
               ),
               tags$li("If desired, clicking \"Copy model to design input\" will transfer the resulting model to the design tab to find the optimal design.")
             ),
             
             sidebarLayout(
               sidebarPanel(
                 "Options",
                 numericInput("fp_bound", "Upper bound", 10, 1, NA, 1),
                 selectInput("fpdegree", "Degree", c(2, 3), selected = 2)
                 
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
             tags$p(
               "This app allows the user to find locally D-optimal designs for fractional polynomials using a variety of metaheuristic optimization algorithms.
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
               tags$li("Click the \"Find\" button to find the optimal design given the inputs.
                       The algorithm should take 10-20 seconds to find the design for default algorithm options.
                       Design points and weights are displayed rounded to 3 decimal places")
               
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
                 actionButton("find", "Find"),
                 waiter::use_waiter(),
                 verbatimTextOutput("design_out")
               )
             )
    )
  )
)