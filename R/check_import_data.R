# check_import_data.R
# make sure user submitted data has the right formatting
# returns the data if correct and errors out if incorrect
# could improve by handling the errors gracefully
# bound: upper bound from interface
check_import_data = function(import_data) {
  d = import_data
  if (ncol(d) != 2)
    stop("Data does not have exactly 2 columns.")
  if (nrow(d) == 0)
    stop("Data is empty")
  if (!identical(colnames(d), c("y", "x"))) {
    warning("Data columns were not named y and x.\nChanging automatically.")
    colnames(d) = c("y", "x")
  }
  if (sum(d$y > 1 | d$y < 0) > 0)
    stop("y column has values outside [0, 1].")
  if (sum(d$x < 0.1) > 0)
    stop("x has values that are smaller than 0.1.")
  
  return(d)
}