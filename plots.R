plot_heatmap <- function(param_matrix){
  
  #check if the matrix is symetric
  if (NROW(param_matrix)==NCOL(param_matrix)){
    param_matrix <- param_matrix-diag(diag(param_matrix))
    param_matrix[upper.tri(param_matrix)]<-NA
  }
  library(reshape2)
  melted_cormat <- melt(param_matrix, na.rm = TRUE)
  # Heatmap
  library(ggplot2)
  ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", space = "Lab") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))
  
}