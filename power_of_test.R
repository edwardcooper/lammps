## Use this script to calculate the power of different stationary test. 


# Generate a matrix of 100*1000 with a mean of 0, standard deviation of 1. 
# I will vary the size of the matrix to see how they affect the power of statistical test. 

library(foreach)
num_row=100
num_column=1000
normal_matrix=foreach(i=1:num_row,.combine=rbind)%do%{
  random_number=rnorm(n=num_column,mean = 0,sd=1)
  return(random_number)
}