library(dplyr)
library(caret)
library(Matrix)

calculate_penalization <- function(matrix_beta,matrix_rho,matrix_phi,lambda,n,p,q,levels_per_variable,t)
{
  # Calculate the penalization term
  # matrix_beta: the matrix of the beta parameters (p,p)
  # matrix_rho: the matrix of the rho parameters (p,q)
  # matrix_phi: the matrix of the phi parameters (q,q)
  # lambda: the lambda parameter
  # n: the number of records
  # p: the number of continuous variables
  # q: the number of discrete variables
  # t: used for thresholding if it is not null
  # return: the penalization term
  
  # Get the cumsum of the levels_per_variable
  cumsum_levels <- cumsum(levels_per_variable)
  
  #Lambda is hardcoded
  lambda <- 5 * sqrt((log(p+q))/n)
  if (!missing(t)){
    lambda <- lambda * t
  }
  
  
  # Do the thresholding if t exists ( Proximal gradient)
  if (!missing(t)){
    beta_weights_aux <- matrix(1,p,p)
    beta_weights <- matrix(pmax(0,beta_weights_aux - t/abs(matrix_beta)),p,p)
    matrix_beta <- matrix_beta * beta_weights
  }
  # Calculate the norm
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
      
      # calculate weights if t exists
      if(!missing(t)){
        weights_aux_matrix_rho <- matrix(1,1,upper_bound - lower_bound + 1)
        weights_aux_matrix_rho <- matrix(pmax(0,weights_aux_matrix_rho - t/abs(aux_matrix_rho)),
                                         1,
                                         upper_bound - lower_bound +1)
        aux_matrix_rho <-aux_matrix_rho * weights_aux_matrix_rho
      }
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
        
        # Apply weights if we need to threshold
        if(!missing(t)){
          weights_aux_matrix_phi <- matrix(1,upper_bound_j - lower_bound_j +1 ,upper_bound_k - lower_bound_k+1)
          weights_aux_matrix_phi <- matrix(pmax(0,weights_aux_matrix_phi - t/abs(aux_matrix_phi)),
                                           upper_bound_j - lower_bound_j +1,
                                           upper_bound_k - lower_bound_k +1) 
          aux_matrix_phi <- aux_matrix_phi * weights_aux_matrix_phi
        }
        
        norm_phi <- norm_phi + norm(aux_matrix_phi, type="F")
      }
      
      
    }
  }
  
  penalization = lambda * (norm_phi+norm_rho+norm_beta)
  return(penalization)
}