rm(list=ls())
library(dplyr)
library(caret)
library(Matrix)
library(ggplot2)


source("matrix_utils.R")
source("preprocessing_utils.R")
source("optimization.R")
source("plots.R")

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

# Calculate the weights as in the paperf
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

lambdas= c(1e-6,
           1e-5,
           1e-4,
           1e-3,
           1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,
           1e-1,5e-1)
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
                   conv=1e-5
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

non_zero_params <- apply(param_df, 2, function(x) sum(x!=0))
df <- data.frame(LogLambda = log10(lambdas),objective= unlist(objectives),non_zero_params=non_zero_params)

ggplot(df, aes(x=LogLambda)) +
  geom_line( aes(y=objective, color="objective") ) + 
  geom_point(aes(y=objective,color="objective"))+
  geom_line( aes(y=non_zero_params, color="non_zero_params")) +
  geom_point(aes(y=non_zero_params,color="non_zero_params"))
  scale_y_continuous(
    name = "Objective",
    sec.axis = sec_axis(~.*1, name="Number of non-zero parameters")
  ) +
  labs(x = "LogLambda", y = "Objective", color = "")+
  scale_color_manual(values = c("red", "blue")) +
  theme(legend.position = "top")
  
# plot the 4th row of the data frame using ggplot
df2 <- data.frame(lambda = log10(lambdas), param = t(as.vector(param_df[4,])))
ggplot(df2, aes(x = lambda, y = X4))  +
   geom_point() +
   geom_line() +
   xlab("LogLambda") +
   ylab("beta_12")





recommended_lambda <- sqrt(log(total_continuous+total_discrete)/total_rows)


param_vec<- convert_params_to_vector(vector_alpha = vector_alpha,
                                     matrix_beta = matrix_beta,
                                     matrix_rho = matrix_rho,
                                     matrix_phi = matrix_phi,
                                     levels_per_variable = levels_per_variable,
                                     p = total_continuous,
                                     q = total_discrete)

results<- proxGD(X = X,Y = Y,
                 initial_param_vector =param_vec,
                 lambda = recommended_lambda,
                 levels_per_variable = levels_per_variable,
                 total_levels = total_levels,
                 n = total_rows,
                 p = total_continuous,
                 q = total_discrete,
                 cont_weights=cont_weights,
                 cat_weights= cat_weights,
                 iter = 5000,
                 gamma=0.0005,
                 conv=1e-5
)

param_list <- convert_vector_to_params(param_vector = results$optimal_params ,
                                       levels_per_variable = levels_per_variable,
                                       p = total_continuous,
                                       q = total_discrete)

alpha <- param_list[[1]]
beta <- param_list[[2]]
rho <- param_list[[3]]
phi <- param_list[[4]]
# change columnames on the matrix rho using the names of Y
colnames(rho) <- colnames(Y)
rownames(rho) <- colnames(X)
# change columnames on the matrix phi using the names of Y
colnames(phi) <- colnames(Y)
rownames(phi) <- colnames(Y)
# change columnames on the matrix beta using the names of X
colnames(beta) <- colnames(X)
rownames(beta) <- colnames(X)


plot_heatmap(beta)
plot_heatmap(rho)
plot_heatmap(phi)
