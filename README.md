# ODpoly

The goal of this software is to provide an easy interface to find optimal designs for logistic regression models with (fractional) polynomial linear predictors. Once installed, call `ODpolyApp()` to launch the R Shiny app. There is also the `ODpoly()` function which accomplishes the  same thing without the user interface. Background information about fractional polynomials and optimal design is included under the "Background" tab in the Shiny app. Usage examples of the R function interface can be found in `ODpoly_example.R`. 

## Try it
You can try a web-based version of the software  [here](https://willgertschapps.shinyapps.io/ODpoly/). For the best performance, I would recommend installing to your local machine. 

## Installation
ODpoly is organized as an R package so the installation can be done directly using devtools.
```
# install.packages("devtools")
devtools::install_github("willgertsch/ODpoly")
```

## List of current features
- Find locally D-optimal designs for given coefficient and powers.
- Supports degree 2 and 3 fractional polynomials.
- Fractional polynomial mini-app that allows the user to fit a fractional polynomial by specifying probability of response data on a graph.
- Uses metaheuristic optimization algorithms from the `metaheuristicOpt` package. The best algorithms for this design problem are selected for use in the Shiny app.
- Automatic merging of identical and 0 weight design points.
