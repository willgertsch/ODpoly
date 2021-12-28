library(shiny)
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
metaheuristics = c("PSO", "ALO", "GWO", "DA", "FFA", "GA", "GOA", "HS", "MFO",
                   "SCA", "WOA", "CLONALG", "DE", "SFL", "CSO", "ABC", "KH", 
                   "CS", "BA", "GBS", "BHO")
ui <- fluidPage(
  titlePanel(
    "Optimal Designs for Fractional Polynomials"
  ),
  sidebarLayout(
    sidebarPanel(
      # inputs
      "Fractional powers:",
      selectInput("p1", "Power 1", frac.powers, selected = 2),
      selectInput("p2", "Power 2", frac.powers, selected = -2),
      "Coefficients:",
      numericInput("b0", "Beta0", 2, -Inf, Inf, 0.01), # bad idea to use inf? probably
      numericInput("b1", "Beta1", 1, -Inf, Inf, 0.01),
      numericInput("b2", "Beta2", -4, -Inf, Inf, 0.01),
      "Algorithm options:",
      selectInput("alg", "Algorithms", metaheuristics, selected = "DE"),
      numericInput("iter", "Iterations", 1000, 1, 10e7, 1),
      numericInput("swarm", "Swarm size", 100, 1, 10e5, 1),
      "Design options:",
      numericInput("pts", "Max design points", 4, 1, 10, 1),
      numericInput("bound", "Upper bound", 10, 1, 10, 1)
    ),
    mainPanel(
      #radioButtons("color", "Pick Color", c("Pink", "Green", "Blue")),
      #selectInput("shape", "Select Shape:", c("Circle", "Triangle")),
      plotOutput("plot1", click = "plot_click"),
      actionButton("fit", "Fit"),
      actionButton("rem_point", "Remove Last Point"),
      actionButton("clear", "Clear all"),
      #verbatimTextOutput("model_out"),
      plotOutput("sens_plot"),
      actionButton("find", "Find optimal design"),
      verbatimTextOutput("design_out")
    )
  )
)