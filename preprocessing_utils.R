library(dplyr)
library(caret)
library(Matrix)

clean_survey_data <- function(data){
  # Drop IDs and region ( has only 1 value)
  drops <- c("X.1","X",'region')
  data <- data[, !(names(data) %in% drops)]
  return(data)
}

get_continuous_discrete <- function(data){
  # keep only the continuous variables
  continuous <- data[, sapply(data, is.numeric)]
  # Convert categorical to factor and save them
  data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], 
                                             as.factor)
  discrete <- data[, sapply(data, is.factor)]
  return(list(continuous, discrete))
}



convert_params_to_vector <- function(vector_alpha,
                         matrix_beta,
                         matrix_rho,
                         matrix_phi,
                         levels_per_variable,
                         p,q){

  total_levels <- sum(levels_per_variable)
  param_sizes <- c(p,p^2,p*total_levels,total_levels^2)
  total_size <-  sum(param_sizes)
  cumsum_sizes <- cumsum(param_sizes)
  
  param_vector=rep(NA,total_size)
  # alpha
  param_vector[1:cumsum_sizes[1]]=vector_alpha
  #beta
  param_vector[(cumsum_sizes[1]+1):cumsum_sizes[2]]<- as.vector(matrix_beta)
  #rho
  param_vector[(cumsum_sizes[2]+1):cumsum_sizes[3]] <- as.vector(matrix_rho)
  #phi
  param_vector[(cumsum_sizes[3]+1):cumsum_sizes[4]] <- as.vector(matrix_phi)
  
  return(param_vector)
  
}


convert_vector_to_params<- function(param_vector,
                                    levels_per_variable,
                                    p,q){
  total_levels <- sum(levels_per_variable)
  param_sizes <- c(p,p^2,p*total_levels,total_levels^2)
  total_size <-  sum(param_sizes)
  cumsum_sizes <- cumsum(param_sizes)
  
  vector_alpha <- param_vector[1:cumsum_sizes[1]]
  matrix_beta <- matrix(param_vector[(cumsum_sizes[1]+1):cumsum_sizes[2]],
                           p,p)
  matrix_rho <- matrix(param_vector[(cumsum_sizes[2]+1):cumsum_sizes[3]],
                       p,total_levels)
  
  matrix_phi <- matrix(param_vector[(cumsum_sizes[3]+1):cumsum_sizes[4]],
                       total_levels,total_levels)
  
  return(list(vector_alpha,matrix_beta,matrix_rho,matrix_phi))
  
}