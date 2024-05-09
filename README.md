# Mixed Graphical Model

Code used for the Modelling under dependence Seminar, based on the work of [Lee and Hastie](https://hastie.su.domains/Papers/structmgm.pdf). Just install the libraries used and source main.R in the mgm-src folder

 - main.R : Main script 
 - optimization.R : Proximal gradient procedure to minimize the objective
 - pseudo_likelihood.R : Calculation of the pseudolikelihood for the continuous and categorical variables
 - penalization.R: functions to Calculate penalty and get the proximal operator
 - gradient.R : get the gradient of the smooth f ( the pseudolikelihood)
 - matrix_utils.R :  common matrix functions used 
 - preprocessing_utils.R : common functions to process data 
 - plots.R : Plot generating functions.


## Documents

- [Presentation](Presentation.pdf)
- [Report](Report.pdf)