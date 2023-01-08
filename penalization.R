library(dplyr)
library(caret)
library(Matrix)


source("preprocessing_utils.R")

calculate_proximal <- function(param_vector, scaled_lambda,p,q,levels_per_variable){
  
  # Get parameters from vector
  param_list <- convert_vector_to_params(param_vector = param_vector,
                                         levels_per_variable = levels_per_variable,
                                         p = p,
                                         q = q)
  vector_alpha<- param_list[[1]]
  matrix_beta<- param_list[[2]]
  matrix_rho<- param_list[[3]]
  matrix_phi  <- param_list[[4]]
  
  cumsum_levels <- cumsum(levels_per_variable)
  
  
  # Do the thresholding for beta
  beta_weights_aux <- matrix(1,p,p)
  beta_weights <- matrix(pmax(0,beta_weights_aux - scaled_lambda/abs(matrix_beta)),p,p)
  matrix_beta <- matrix_beta * beta_weights
  
  # Scale rho matrix
  # For all continuous
  for (i in 1:p){
    # For all discrete
    for (j in 1:q){
      levels_in_this_varible = levels_per_variable[j]
      cummulated_levels_in_this_variable <- cumsum_levels[j]
      lower_bound <- cummulated_levels_in_this_variable-levels_in_this_varible+1
      upper_bound <- cummulated_levels_in_this_variable
      aux_matrix_rho <- matrix_rho[i,lower_bound:upper_bound]
      
      # calculate weights if t exists
      weights_aux_matrix_rho <- matrix(1,1,upper_bound - lower_bound + 1)
      weights_aux_matrix_rho <- matrix(pmax(0,weights_aux_matrix_rho - scaled_lambda/abs(aux_matrix_rho)),
                                       1,
                                       upper_bound - lower_bound +1
                                       )
      aux_matrix_rho <-aux_matrix_rho * weights_aux_matrix_rho
      
      matrix_rho[i,lower_bound:upper_bound] <- aux_matrix_rho
      }
    }
  
    # Scale phi matrix
    for (j in 1:q){
      levels_in_j = levels_per_variable[j]
      cummulated_levels_in_j <- cumsum_levels[j]
      #For all discrete greater than j and less than q
      for (k in (j+1):q){
        if (k<= q) {
          levels_in_k = levels_per_variable[k]
          cummulated_levels_in_k <- cumsum_levels[k]
          lower_bound_j <- cummulated_levels_in_j-levels_in_j+1
          lower_bound_k <- cummulated_levels_in_k-levels_in_k+1
          upper_bound_j <- cummulated_levels_in_j
          upper_bound_k <- cummulated_levels_in_k
          
          aux_matrix_phi <- matrix_phi[lower_bound_j:upper_bound_j,lower_bound_k:upper_bound_k]
          
          # Apply threshold
          weights_aux_matrix_phi <- matrix(1,upper_bound_j - lower_bound_j +1 ,upper_bound_k - lower_bound_k+1)
          weights_aux_matrix_phi <- matrix(pmax(0,weights_aux_matrix_phi - scaled_lambda/abs(aux_matrix_phi)),
                                           upper_bound_j - lower_bound_j +1,
                                           upper_bound_k - lower_bound_k +1) 
          aux_matrix_phi <- aux_matrix_phi * weights_aux_matrix_phi
          
          matrix_phi[lower_bound_j:upper_bound_j,lower_bound_k:upper_bound_k] <- aux_matrix_phi
          }
        }
    }
  
  # Create vector from gradients
  threshold_vector<- convert_params_to_vector(vector_alpha = vector_alpha,
                                          matrix_beta = matrix_beta,
                                          matrix_rho = matrix_rho,
                                          matrix_phi = matrix_phi,
                                          levels_per_variable = levels_per_variable,
                                          p = p,
                                          q = q) 
  
  return(threshold_vector)
  
  }
  
  
  
  
  

calculate_penalization <- function(param_vector,lambda,n,p,q,levels_per_variable)
{
  # Calculate the penalization term
  # param_vector: a vector that can be converted to all the parameters
  # lambda: the lambda parameter
  # n: the number of records
  # p: the number of continuous variables
  # q: the number of discrete variables
  # t: used for thresholding if it is not null
  # return: the penalization term
  
  param_list <- convert_vector_to_params(param_vector = param_vector,
                                         levels_per_variable = levels_per_variable,
                                         p = p,
                                         q = q)
  vector_alpha<- param_list[[1]]
  matrix_beta<- param_list[[2]]
  matrix_rho<- param_list[[3]]
  matrix_phi  <- param_list[[4]]
  
  # Get the cumsum of the levels_per_variable
  cumsum_levels <- cumsum(levels_per_variable)
  
  # Calculate the norm of beta
  norm_beta <- sum(abs(matrix_beta))
  
  # Calculate the penalization term for rho matrix
  norm_rho <- 0
  # For all continuous
  for (i in 1:p){
    # For all discrete
    for (j in 1:q){
      levels_in_this_varible = levels_per_variable[j]
      cummulated_levels_in_this_variable <- cumsum_levels[j]
      lower_bound <- cummulated_levels_in_this_variable-levels_in_this_varible+1
      upper_bound <- cummulated_levels_in_this_variable
      aux_matrix_rho <- matrix_rho[i,lower_bound:upper_bound]
      
      norm_rho <- norm_rho + norm(aux_matrix_rho, type="2")
    }
  }
  
  
  # Calculate the penalization term for the phi matrix
  norm_phi<-0
  # For all discrete
  for (j in 1:q){
    levels_in_j = levels_per_variable[j]
    cummulated_levels_in_j <- cumsum_levels[j]
    #For all discrete greater than j and less than q
    for (k in (j+1):q){
      if (k<= q) {
        levels_in_k = levels_per_variable[k]
        cummulated_levels_in_k <- cumsum_levels[k]
        lower_bound_j <- cummulated_levels_in_j-levels_in_j+1
        lower_bound_k <- cummulated_levels_in_k-levels_in_k+1
        upper_bound_j <- cummulated_levels_in_j
        upper_bound_k <- cummulated_levels_in_k
        
        aux_matrix_phi <- matrix_phi[lower_bound_j:upper_bound_j,lower_bound_k:upper_bound_k]
        norm_phi <- norm_phi + norm(aux_matrix_phi, type="F")
      }
      
      
    }
  }
  
  
  
  penalization = lambda * (norm_phi+norm_rho+norm_beta)
  return(penalization)
}