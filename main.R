rm(list=ls())
library(dplyr)
library(caret)
library(Matrix)


source("matrix_utils.R")
source("preprocessing_utils.R")
source("optimization.R")

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

# Calculate the weights as in the paper
cont_weights <- matrix(sqrt(sapply(X, var)),1,NCOL(X))
relative_frequencies <- sapply(Y, function(x) table(x)/length(x))
cat_weights <- matrix(0,1,NCOL(Y))
for (i in 1:length(relative_frequencies)){
  aux<- unlist(relative_frequencies[i])
  cat_weights[i]<-  sum(sqrt(aux * (1 - aux)))
}


# Number of continouous variables
total_continuous <- NCOL(X)
# Number of discrete variables
total_discrete <- NCOL(Y)
# Rows
total_rows <- NROW(X)



# Convert dataframe to matrix
X<- data.matrix(X)
#do one-hot encoding to the discrete variables\
Y<- dummyVars("~ .", data = Y) %>% predict(Y)
Y<- data.matrix(Y)

lambdas= c(1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5)
objectives= list()
param_df = NA

matrix_beta <- create_simmmetryc_matrix(total_continuous,total_continuous)
diag(matrix_beta) <- 1
vector_alpha <- create_simmmetryc_matrix(1,total_continuous) 
matrix_rho <- create_simmmetryc_matrix(total_continuous,total_levels)
matrix_phi <- create_simmmetryc_matrix(total_levels,total_levels)

for (lambda in lambdas){
  print(lambda)
  param_vec<- convert_params_to_vector(vector_alpha = vector_alpha,
                                       matrix_beta = matrix_beta,
                                       matrix_rho = matrix_rho,
                                       matrix_phi = matrix_phi,
                                       levels_per_variable = levels_per_variable,
                                       p = total_continuous,
                                       q = total_discrete)
  
  results<- proxGD(X = X,Y = Y,
                   initial_param_vector =param_vec,
                   lambda = lambda,
                   levels_per_variable = levels_per_variable,
                   total_levels = total_levels,
                   n = total_rows,
                   p = total_continuous,
                   q = total_discrete,
                   cont_weights=cont_weights,
                   cat_weights= cat_weights,
                   iter = 5000,
                   gamma=0.0005,
                   conv=1e-4
  )
  param_list <- convert_vector_to_params(param_vector = results$optimal_params ,
                                         levels_per_variable = levels_per_variable,
                                         p = total_continuous,
                                         q = total_discrete)
  
  if (is.na(param_df)){
    param_df = data.frame(results$optimal_params)
    }
  else{
    param_df[,ncol(param_df)+1]<- results$optimal_params
  }
  objectives <- append(objectives, results$objective[length(results$objective)])
}

plot(log(lambdas),objectives)
