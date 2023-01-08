library(dplyr)
library(caret)
library(Matrix)


force_simmetry <- function(input_matrix){
  ind <- lower.tri(input_matrix)
  input_matrix[ind] <- t(input_matrix)[ind]
  return (input_matrix)
}


create_simmmetryc_matrix<- function(nrows,ncols){
  mat <- matrix(rnorm(ncols*nrows,mean=0,sd=0.1),
                nrows,
                ncols)
  
  mat <- force_simmetry(mat)
  return(mat)
}