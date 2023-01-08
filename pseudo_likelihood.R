library(dplyr)
library(caret)
library(Matrix)

source("matrix_utils.R")

calculate_discrete_pseudo_likelihood <- function(X,Y,
                                                 matrix_rho,
                                                 matrix_phi,
                                                 levels_per_variable,
                                                 total_levels,
                                                 n,
                                                 p,
                                                 q){
  # Calculate the discrete likelihood
  # Y: the discrete variables with one-hot encoding applied (n,total_levels)
  # matrix_phi: the matrix of the phi parameters (q,q)
  # levels_per_variable: the number of levels of each variable (q)
  # total_levels: the total number of levels (q)
  # return: the discrete pseudo-likelihood
  
  cumsum_levels <- cumsum(levels_per_variable)
  discrete_pseudo_likelihood <-0
  
  # get the diagonal of phi matrix and substract it from phi
  diag_phi <- diag(matrix_phi)
  matrix_phi_no_diag <- matrix_phi - diag(diag(matrix_phi))
  # create a vectors of ones to do the sum as a matrix multiplication
  e <- rep(1, n)
  
  # calculate the inside part of the softmax for each variable for each record
  matrix_w <- X %*% matrix_rho + Y %*% matrix_phi_no_diag + e %*% matrix(diag_phi,1)
  
  for (i in 1:q){
    # get the number of levels of the variable
    levels_in_this_varible = levels_per_variable[i]
    cummulated_levels_in_this_variable = cumsum_levels[i]
    lower_bound = cummulated_levels_in_this_variable-levels_in_this_varible+1
    upper_bound = cumsum_levels[i]
    
    # get the inside part of the softmax for the discrete variable
    aux_matrix_w <- matrix_w[,lower_bound:upper_bound]
    aux_matrix_Y <- Y[,lower_bound:upper_bound]
    # calculate the softmax denominator for each row ( using all levels)
    # TODO: Change to use the correct column
    aux_denominator <- apply(aux_matrix_w, 1, function(x) sum(exp(x)))
    
    matrix_numerators= exp(aux_matrix_w)*aux_matrix_Y
    aux_numerator <- apply(matrix_numerators, 1, function(x) sum(x))
    
    discrete_pseudo_likelihood <- discrete_pseudo_likelihood - sum(log(aux_numerator/aux_denominator))
  }
  
  return(discrete_pseudo_likelihood)
  
  
}

calculate_continuous_pseudo_likelihood <- function(X,Y,
                                                   matrix_beta,
                                                   vector_alpha,
                                                   matrix_rho,
                                                   n,
                                                   p,
                                                   q){
  # Calculate the continuous pseudo likelihood
  # X: the continuous variables (n,p)
  # Y: the discrete variables with one-hot encoding applied (n,total_levels)
  # matrix_beta: the matrix of the beta parameters (p,p)
  # vector_alpha: the vector of the alpha parameters (p)
  # return: the continuous likelihood
  
  # get the diagonal of beta matrix to calculate the first term ( beta_ss)
  diag_beta <- diag(matrix_beta)
  # create a vectors of ones to do the sum as a matrix multiplication
  e <- rep(1, n)
  
  # First term of the likelihood
  continuous_pseudo_likelihood <- - n / 2 * sum(log(diag_beta))
  
  # substract the diagonal to beta
  matrix_beta_no_diag <- matrix_beta - diag(diag(matrix_beta))
  matrix_beta_no_diag<- force_simmetry(matrix_beta_no_diag)
  X_no_diag <- X %*% matrix_beta_no_diag %*% diag(1/diag_beta)
  Y_rho <- Y %*% t(matrix_rho) %*%  diag(1/diag_beta)
  # Alpha Component
  alpha_component <- e %*% vector_alpha %*% diag(1/diag_beta)
  
  # Calculate the second term of the continuous pseudo likelihood
  continuous_pseudo_likelihood <- continuous_pseudo_likelihood +
    .5 * norm((X - alpha_component + X_no_diag - Y_rho) %*% diag(sqrt(diag_beta)),
              type = 'F')^2
  
  # Issue here, the matlab code has different sigs and the alpha are not divided by diag_beta
  
  return(continuous_pseudo_likelihood)
}

calculate_pseudo_likelihood <- function(X,Y,
                                        matrix_beta,
                                        vector_alpha,
                                        matrix_rho,
                                        matrix_phi,
                                        levels_per_variable,total_levels,
                                        n,p,q){
  # Calculate the pseudo likelihood
  # X: the continuous variables (n,p)
  # Y: the discrete variables (n,q)
  # matrix_beta: the matrix of the beta parameters (p,p)
  # vector_alpha: the vector of the alpha parameters (p)
  # matrix_rho: the matrix of the rho parameters (p,q)
  # matrix_phi: the matrix of the phi parameters (q,q)
  # return: the pseudo likelihood
  # Number of continouous variables
  
  # Force simmetry 
  matrix_beta <- force_simmetry(matrix_beta)
  matrix_phi <- force_simmetry(matrix_phi)
  

  
  continuous_pseudo_likelihood <- calculate_continuous_pseudo_likelihood(X,Y,
                                                                         matrix_beta,
                                                                         vector_alpha,
                                                                         matrix_rho,
                                                                         n = n,
                                                                         p = p,
                                                                         q = q)
  
  # Calculate the discrete likelihood
  discrete_pseudo_likelihood <- calculate_discrete_pseudo_likelihood(X,Y,
                                                                     matrix_rho = matrix_rho,
                                                                     matrix_phi = matrix_phi,
                                                                     levels_per_variable = levels_per_variable,
                                                                     total_levels = total_levels,
                                                                     n = n,
                                                                     p = p,
                                                                     q = q
                                                                     )
  
  return((continuous_pseudo_likelihood+discrete_pseudo_likelihood)/n)
  
}




