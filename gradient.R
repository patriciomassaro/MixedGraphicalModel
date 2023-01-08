library(dplyr)
library(caret)
library(Matrix)

source ("matrix_utils.R")
source("preprocessing_utils.R")

calculate_gradient_beta<- function(matrix_beta_no_diag,diag_beta,
                                   matrix_rho,
                                   vector_alpha,
                                   X,
                                   X_beta,
                                   Y_rho,
                                   res,
                                   n,p,q){
    #compute the gradient
    gradient_beta = t(X) %*% res
    # zero the diagonal and force simmetry
    gradient_beta_diag <- diag(gradient_beta)
    gradient_beta_no_diag <- force_simmetry (gradient_beta - diag(gradient_beta_diag))

    # add the diagonal component
    gradient_beta_diag <- rep(0, p)
    for (i in 1:p){
        gradient_beta_diag[i] <- - ( n / (2*diag_beta[i]) ) + 1/2*norm(res[,i],type = '2')^2 - res[,i] %*% (X_beta[,i] + Y_rho[,i])
    }

    # add the diagonal component to the gradient
    gradient_beta <- gradient_beta_no_diag + diag(gradient_beta_diag)

    return(gradient_beta)


}



calculate_gradient_alpha <- function(matrix_beta,res){
  grad_alpha <- diag(matrix_beta) * t(matrix(colSums(res)))
  return (grad_alpha)
  
}

calculate_gradient_rho_and_phi <- function(X,Y,
                                           res,
                                           matrix_w,
                                           levels_per_variable,cumsum_levels,q){
  grad_rho_first_term <- t(t(Y) %*%(res))
  
  for(i in 1:q){
    levels_in_this_varible = levels_per_variable[i]
    cummulated_levels_in_this_variable = cumsum_levels[i]
    lower_bound = cummulated_levels_in_this_variable-levels_in_this_varible+1
    upper_bound = cumsum_levels[i]
    
    aux_matrix_w <- matrix_w[,lower_bound:upper_bound]
    
    aux_denominator <- apply(aux_matrix_w, 1, function(x) sum(exp(x)))
    
    aux_matrix_w = diag(1/aux_denominator) %*% exp(aux_matrix_w) - 1
    
    matrix_w[,lower_bound:upper_bound] <- aux_matrix_w
  }
  
  
  
  grad_rho_second_term <- t(X) %*% matrix_w
  
  grad_rho = grad_rho_first_term + grad_rho_second_term
  
  
  grad_phi <- t(Y) %*% matrix_w
  grad_phi <- force_simmetry(grad_phi)
  diag(grad_phi) <- colSums(matrix_w)
  
  return(list(grad_rho,grad_phi))
}

calculate_gradient <- function(param_vector,
                               X,Y,levels_per_variable,
                               p,q,n){
# Calculates the gradient of the pseudo-likelihood function ( our smooth F)
#  param_vector: p + p*p + p*q + q*q vector
# X: n x p matrix
# Y: n x total_levels matrix
# levels_per_variable: 1 x q matrix
# n: number of rows
# p: number of continuous variables
# q: number of discrete variables
  
  
  param_list <- convert_vector_to_params(param_vector = param_vector,
                                         levels_per_variable = levels_per_variable,
                                         p = p,
                                         q = q)
  vector_alpha<- matrix(param_list[[1]],1,p)
  matrix_beta<- param_list[[2]]
  matrix_rho<- param_list[[3]]
  matrix_phi  <- param_list[[4]]
  
  total_levels = sum(levels_per_variable)
  cumsum_levels = cumsum(levels_per_variable)

  # substract the diagonal of the beta matrix and force symmetry
  diag_beta <- diag(matrix_beta)
  matrix_beta_no_diag <- force_simmetry (matrix_beta - diag(diag_beta))

  # substract the diagonal of the phi matrix and force symmetry
  diag_phi <- diag(matrix_phi)
  matrix_phi_no_diag <- force_simmetry (matrix_phi - diag(diag_phi))

  # Prepare the quantities used by the gradients only once
  X_beta = X %*% matrix_beta_no_diag %*%  diag(1/diag_beta)
  Y_rho = Y %*% t(matrix_rho) %*% diag(1/diag_beta)
  e <- rep(1, n)
  alpha_component = e %*% vector_alpha 
  res <- X_beta - X + alpha_component + Y_rho
  matrix_w <- X %*% matrix_rho + Y %*% matrix_phi_no_diag + e %*% matrix(diag_phi,1)
  
  
  # Calculate each gradient
  grad_beta <- calculate_gradient_beta(matrix_beta_no_diag = matrix_beta_no_diag,
                                       diag_beta = diag_beta,
                                       matrix_rho =matrix_rho,
                                       vector_alpha = vector_alpha,
                                       X = X,
                                       X_beta = X_beta,
                                       Y_rho = Y_rho,
                                       res = res,
                                       n=n,p=p,q=q
                                       )/n

  grad_alpha <- calculate_gradient_alpha(matrix_beta = matrix_beta,
                                         res=res)/n
  
  
  aux <- calculate_gradient_rho_and_phi(X = X,Y = Y,res = res,
                                        matrix_w = matrix_w,
                                        levels_per_variable = levels_per_variable,
                                        cumsum_levels = cumsum_levels,
                                        q = q
                                        )
  # Expand the list with rho and phi
  grad_rho <- aux[[1]]/n
  grad_phi <- aux[[2]]/n
  
  
  # Create vector from gradients
  gradient_vec<- convert_params_to_vector(vector_alpha = grad_alpha,
                                       matrix_beta = grad_beta,
                                       matrix_rho = grad_rho,
                                       matrix_phi = grad_phi,
                                       levels_per_variable = levels_per_variable,
                                       p = p,
                                       q = q)  


  return(gradient_vec)
}

