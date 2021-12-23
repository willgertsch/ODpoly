library(shiny)

# globals
frac.powers = c(-2, -1, -1/2, 0, 1/2, 1, 2, 3)
xs = seq(0, 10, 0.01)

ODpolyApp <- function(...) {
  ui <- fluidPage(
    titlePanel(
      "Optimal Designs for Fractional Polynomials"
    ),
    sidebarLayout(
      sidebarPanel(
        # inputs
        "Fractional powers:",
        selectInput("p1", "Power 1", frac.powers, selected = 1),
        selectInput("p2", "Power 2", frac.powers, selected = 2),
        "Coefficients:",
        numericInput("b0", "Beta0", 1, -Inf, Inf, 0.01), # bad idea to use inf? probably
        numericInput("b1", "Beta1", 1, -Inf, Inf, 0.01),
        numericInput("b2", "Beta2", 1, -Inf, Inf, 0.01)
      ),
      mainPanel(
        #radioButtons("color", "Pick Color", c("Pink", "Green", "Blue")),
        #selectInput("shape", "Select Shape:", c("Circle", "Triangle")),
        actionButton("rem_point", "Remove Last Point"),
        actionButton("clear", "Clear all"),
        plotOutput("plot1", click = "plot_click"),
        actionButton("fit", "Find fractional poly"),
        #verbatimTextOutput("model_out")
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
        lims(x = c(0, 10), y = c(0, 1)) +
        theme_bw() + 
        # include so that colors don't change as more color/shape chosen
        #scale_color_discrete(drop = FALSE) +
        #scale_shape_discrete(drop = FALSE) +
        labs(y = "probability of response", x = "dose")
      
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
      add_row <- data.frame(x = input$plot_click$x,
                            y = input$plot_click$y,
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
  }
  
  
  shinyApp(ui, server, ...)
}