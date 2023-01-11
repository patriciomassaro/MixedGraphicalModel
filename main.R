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

# Initiate parameters
matrix_beta <- create_simmmetryc_matrix(total_continuous,total_continuous)
diag(matrix_beta) <- 1
vector_alpha <- create_simmmetryc_matrix(1,total_continuous) # CONT
matrix_rho <- create_simmmetryc_matrix(total_continuous,total_levels)
matrix_phi <- create_simmmetryc_matrix(total_levels,total_levels)
lambda <-  0.01

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
                 iter = 15
                 )

param_list <- convert_vector_to_params(param_vector = results$optimal_params ,
                                       levels_per_variable = levels_per_variable,
                                       p = total_continuous,
                                       q = total_discrete)
final_vector_alpha<- param_list[[1]]
final_matrix_beta<- param_list[[2]]
final_matrix_rho<- param_list[[3]]
final_matrix_phi  <- param_list[[4]]


print(results$objective)
