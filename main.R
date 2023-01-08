library(dplyr)
library(caret)
library(Matrix)


source("matrix_utils.R")
source("preprocessing_utils.R")
source("pseudo_likelihood.R")
source("penalization.R")
source("gradient.R")

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


# Convert dataframe to matrix
X<- data.matrix(X)
#do one-hot encoding to the discrete variables\
Y<- dummyVars("~ .", data = Y) %>% predict(Y)
Y<- data.matrix(Y)

l <- calculate_pseudo_likelihood(X,Y,
                                 matrix_beta = matrix_beta,
                                 vector_alpha = vector_alpha,
                                 matrix_rho = matrix_rho,
                                 matrix_phi = matrix_phi,
                                 levels_per_variable = levels_per_variable,
                                 total_levels = total_levels,
                                 n= total_rows,
                                 p= total_continuous,
                                 q= total_discrete)
print(l)
l <- l+ calculate_penalization(matrix_beta = matrix_beta,
                               matrix_rho = matrix_rho,
                               matrix_phi = matrix_phi,
                               p = total_continuous,
                               q = total_discrete,
                               n = total_rows,
                               levels_per_variable = levels_per_variable
                               #t=1
                               )
print(l)

grad <- calculate_gradient(matrix_beta = matrix_beta,
                           matrix_rho = matrix_rho,
                           matrix_phi = matrix_phi,
                           vector_alpha = vector_alpha,
                           X,Y,levels_per_variable,
                           p = total_continuous,
                           q = total_discrete,
                           n = total_rows)

grad_beta = grad[[1]]

