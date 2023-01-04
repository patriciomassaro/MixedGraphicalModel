library(dplyr)
library(caret)
library(Matrix)


clean_survey_data <- function(data){
  # Drop IDs and region ( has only 1 value)
  drops <- c("X.1","X",'region')
  data <- data[, !(names(data) %in% drops)]
  return(data)
}

force_simmetry <- function(input_matrix){
  ind <- lower.tri(input_matrix)
  input_matrix[ind] <- t(input_matrix)[ind]
  return (input_matrix)
  
}

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
  # return: the penalization term

  # Get the cumsum of the levels_per_variable
  cumsum_levels <- cumsum(levels_per_variable)
  
  #Lambda is hardcoded
  lambda <- 5 * sqrt((log(p+q))/n) * t

  # Do the thresholding ( Proximal gradient)
  beta_weights_aux <- matrix(1,p,p)
  beta_weights <- matrix(pmax(0,beta_weights_aux - t/abs(matrix_beta)),p,p)
  matrix_beta <- matrix_beta * beta_weights
  
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
      
      # calculate weights
      weights_aux_matrix_rho <- matrix(1,1,upper_bound - lower_bound + 1)
      weights_aux_matrix_rho <- matrix(pmax(0,weights_aux_matrix_rho - t/abs(aux_matrix_rho)),
                               1,
                               upper_bound - lower_bound +1)
      aux_matrix_rho <-aux_matrix_rho * weights_aux_matrix_rho
      
      
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
        # Apply weights
        weights_aux_matrix_phi <- matrix(1,upper_bound_j - lower_bound_j +1 ,upper_bound_k - lower_bound_k+1)
        weights_aux_matrix_phi <- matrix(pmax(0,weights_aux_matrix_phi - t/abs(aux_matrix_phi)),
                                         upper_bound_j - lower_bound_j +1,
                                         upper_bound_k - lower_bound_k +1) 
        aux_matrix_phi <- aux_matrix_phi * weights_aux_matrix_phi
        norm_phi <- norm_phi + norm(aux_matrix_phi, type="F")
      }

      
    }
  }
    
  penalization = lambda * (norm_phi+norm_rho+norm_beta)
    
  return(penalization)
}

create_simmmetryc_matrix<- function(nrows,ncols){
  mat <- matrix(rnorm(ncols*nrows,mean=0,sd=0.1),
                nrows,
                ncols)
  
  mat <- force_simmetry(mat)
  return(mat)
}

calculate_discrete_pseudo_likelihood <- function(X,Y,matrix_rho,matrix_phi,levels_per_variable, total_levels){
  # Calculate the discrete likelihood
  # Y: the discrete variables with one-hot encoding applied (n,total_levels)
  # matrix_phi: the matrix of the phi parameters (q,q)
  # levels_per_variable: the number of levels of each variable (q)
  # total_levels: the total number of levels (q)
  # return: the discrete pseudo-likelihood

  # Number of samples
  n <- NROW(Y)
  q <- length(levels_per_variable)
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
    # calculate the softmax denominator for each row ( using all levels)
    aux_denominator <- apply(aux_matrix_w, 1, function(x) sum(exp(x)))
    for (j in 1:levels_per_variable[i]){
      # calculate the softmax numerator for each row ( one for each level)
      aux_numerator <- exp(aux_matrix_w[,j])
      # calculate the likelihood, summing for every record
      discrete_pseudo_likelihood <- discrete_pseudo_likelihood - sum(log(aux_numerator/aux_denominator))
      
    }
  }
  
  return(discrete_pseudo_likelihood)


}

calculate_continuous_pseudo_likelihood <- function(X,Y, matrix_beta, vector_alpha,matrix_rho){
  # Calculate the continuous pseudo likelihood
  # X: the continuous variables (n,p)
  # Y: the discrete variables with one-hot encoding applied (n,total_levels)
  # matrix_beta: the matrix of the beta parameters (p,p)
  # vector_alpha: the vector of the alpha parameters (p)
  # return: the continuous likelihood


  # Number of records
  n <- NROW(X)
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
  alpha_component <- e %*% vector_alpha 

  # Calculate the second term of the continuous pseudo likelihood
  continuous_pseudo_likelihood <- continuous_pseudo_likelihood +
    .5 * norm((X - alpha_component - X_no_diag - Y_rho) %*% diag(sqrt(diag_beta)),
               type = 'F')^2
  
  return(continuous_pseudo_likelihood)
}

calculate_pseudo_likelihood <- function(X,Y,
                                        matrix_beta,
                                        vector_alpha,
                                        matrix_rho,
                                        matrix_phi,
                                        levels_per_variable,
                                        total_levels){
  # Calculate the pseudo likelihood
  # X: the continuous variables (n,p)
  # Y: the discrete variables (n,q)
  # matrix_beta: the matrix of the beta parameters (p,p)
  # vector_alpha: the vector of the alpha parameters (p)
  # matrix_rho: the matrix of the rho parameters (p,q)
  # matrix_phi: the matrix of the phi parameters (q,q)
  # return: the pseudo likelihood
  # Number of continouous variables
  p <- NCOL(X)
  q <- NCOL(Y)
  n <- NROW(X)
  
  # Force simmetry 
  matrix_beta <- force_simmetry(matrix_beta)
  matrix_rho <- force_simmetry(matrix_rho)
  matrix_phi <- force_simmetry(matrix_phi)
  
  # Convert dataframe to matrix
  X<- data.matrix(X)
  #do one-hot encoding to the discrete variables\
  Y<- dummyVars("~ .", data = Y) %>% predict(Y)
  Y<- data.matrix(Y)

  continuous_pseudo_likelihood <- calculate_continuous_pseudo_likelihood(X,
                                                                  Y,
                                                                  matrix_beta,
                                                                  vector_alpha,
                                                                  matrix_rho)
  
  # Calculate the discrete likelihood
  discrete_pseudo_likelihood <- calculate_discrete_pseudo_likelihood(X,Y,
                                                              matrix_rho,
                                                              matrix_phi,
                                                              levels_per_variable,
                                                              total_levels
                                                              )

  return((continuous_pseudo_likelihood+discrete_pseudo_likelihood)/n)

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


#load the data
data <- read.csv("data/Mid-Atlantic_Wage_Data_974_40.csv")
# convert year column to str
data[,3] <- as.character(data[,3])
data <- clean_survey_data(data)
data <- get_continuous_discrete(data)

# Divide the dataframes into continuous and discrete
X <- data[[1]]
Y <- data[[2]]
# get the numbers of levels of each factor 
levels_per_variable <- sapply(Y, function(x) length(levels(x)))
total_levels <- sum(levels_per_variable)

# Number of continouous variables
total_continuous <- NCOL(X)
# Number of discrete variables
total_discrete <- NCOL(Y)
# Rows
total_rows <- NROW(X)

# Initiate parameters
matrix_beta <- create_simmmetryc_matrix(total_continuous,total_continuous)
matrix_beta <- diag(sign(diag(matrix_beta))) %*% matrix_beta
vector_alpha <- create_simmmetryc_matrix(1,total_continuous) # CONT
matrix_rho <- create_simmmetryc_matrix(total_continuous,total_levels)
matrix_phi <- create_simmmetryc_matrix(total_levels,total_levels)

l <- calculate_pseudo_likelihood(X,Y,matrix_beta,vector_alpha,matrix_rho,matrix_phi,levels_per_variable,total_levels)
print(l)
l <- l+ calculate_penalization(matrix_beta = matrix_beta,
                               matrix_rho = matrix_rho,
                               matrix_phi = matrix_phi,
                               p = total_continuous,
                               q = total_discrete,
                               n = total_rows,
                               levels_per_variable = levels_per_variable,
                               t=1
                               )

print(l)

