library(shiny)
library(waiter)

# globals
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
metaheuristics = c("Particle Swarm Optimization", 
                   "Grey Wolf Optimizer", 
                   "Harmony Search Algorithm", 
                   "Moth Flame Optimizer",
                   "Differential Evolution")

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