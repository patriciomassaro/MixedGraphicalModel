
source("penalization.R")
source("preprocessing_utils.R")
source("pseudo_likelihood.R")
source("gradient.R")

proxGD <- function(X,Y,
                   initial_param_vector,
                   lambda,
                   levels_per_variable,
                   total_levels,
                   n,p,q,
                   gamma=0.0001,
                   iter=5,
                   conv=1e-4){
    # X: n x p matrix
    # Y: n x q matrix
    # initial_param_vector: p + p*p + p*q + q*q vector
    # lambda: scalar (penalization parameter)
    # gamma: scalar (step size)
    # iter: scalar (max number of iterations)
    # conv: scalar (convergence threshold)
  
  
  # Initialize coefficients matrix 
  parameters <- matrix(rep(NA,
                           length(initial_param_vector)*(iter+1)),
                           nrow=length(initial_param_vector)
                       )
  parameters[, 1] <- initial_param_vector

  # Initialize objective funtion
  obj <- rep(NA, iter+1)
  obj[1] <- calculate_pseudo_likelihood(X = X,Y = Y,
                                        param_vec = parameters[,1],
                                        levels_per_variable = levels_per_variable,
                                        total_levels = total_levels,
                                        n = n,
                                        p = p,
                                        q = q
                                        ) 
  obj[1] <- obj[1] + calculate_penalization(param_vector = initial_param_vector,
                                            lambda = lambda,
                                            n = n,p = p,q = q,
                                            levels_per_variable = levels_per_variable
                                            )
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  for (t in 1:iter) {
    gradient <- calculate_gradient(param_vector = parameters[,t],
                                   X = X,Y = Y,
                                   levels_per_variable = levels_per_variable,
                                   p = p,q = q,n = n)
    
    u <- parameters[,t] - gamma * gradient
    parameters[,t+1] <- calculate_proximal(param_vector = u,
                                           scaled_lambda = gamma*lambda,
                                           levels_per_variable = levels_per_variable,
                                           p = p,q = q)
    
    obj[t+1] <- calculate_pseudo_likelihood(X = X,Y = Y,
                                            param_vec =  parameters[,(t+1)],
                                            levels_per_variable = levels_per_variable,
                                            total_levels = total_levels,
                                            n = n,
                                            p = p,
                                            q = q
                                            )
    obj[t+1]<- obj[t+1]+ calculate_penalization(param_vector = parameters[,t+1],
                                                lambda = lambda,
                                                levels_per_variable = levels_per_variable,
                                                n = n,p = p,q = q
                                                )
    
    # Check convergence
    delta <- abs(obj[t+1]-obj[t]) /(abs(obj[t])+conv)
    if (delta < conv) {
      # Remove excess parameters and objectives
      parameters <- parameters[, -((t+2):ncol(parameters))]
      obj <- obj[-((t+2):length(obj))]
      
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t+1))
      break
    }
  }
  
  result <- list("optimal_params"=parameters[,ncol(parameters)],
                 "params_history"=parameters,
                 "objective"=obj,
                 "conv"=message)
  return(result)
}