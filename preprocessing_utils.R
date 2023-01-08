library(dplyr)
library(caret)
library(Matrix)

clean_survey_data <- function(data){
  # Drop IDs and region ( has only 1 value)
  drops <- c("X.1","X",'region')
  data <- data[, !(names(data) %in% drops)]
  return(data)
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