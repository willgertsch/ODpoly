# ODpoly

The goal of this software is to provide an easy interface to find optimal designs for logistic regression models with (fractional) polynomial linear predictors. Once installed, call ODpolyApp() to launch the R Shiny app. There is also the ODpoly() function which accomplishes the  same thing without the user interface. At the time of writing, I haven't written documentation so take a look at ODpoly.R to get an idea of what the arguments should be.

## Try it
You can try a web-based version of the software  [here](https://willgertschapps.shinyapps.io/ODpoly/). For the best performance, I would recommend installing to your local machine.

## Installation
ODpoly is organized as an R package so the installation should be simple.
```
# install.packages("devtools")
devtools::install_github("willgertsch/ODpoly")
```

## List of current features
- Find local D-optimal designs for given coefficient and powers.
- Can click on top plot to generate data for a fractional polynomial model.
- Implements all algorithms from the `metaheuristicOpt` package.
- Automatic merging of identical and 0 weight design points.
