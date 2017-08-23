## Use this script to calculate the power of different stationary test. 


# Generate a matrix of 100*1000 with a mean of 0, standard deviation of 1. 
# I will vary the size of the matrix to see how they affect the power of statistical test. 

# due to the limit on the p value produced the type I error we could choose from is (0.01,0.1)

# write the p_value of statistical test into a function P_value_matrix(). 
P_value_matrix=function(num_test=10,num_sample=1000){
  library(foreach)
 
  p_value_matrix=foreach(i=1:num_test,.combine=rbind)%do%{
    # create vector with a length of num_sample
    # The random vector is a gaussian distribution with mean 0,sd 1. 
    random_number=rnorm(n=num_sample,mean = 0,sd=1)
    
    # do the stationarity test from fractal package, then extract the p value at the end
    test_fractal=fractal::stationarity(random_number)
    
    # do adf test from tseries package
    tseries_adf=tseries::adf.test(random_number)
    # warnings()%>%print()
    
    # do kpss test from tseries package
    tseries_kpss=tseries::kpss.test(random_number,null="Level")
    
    
    # do pp test from tseries package
    tseries_pp= tseries::pp.test(random_number)
    
    #do pp test from the base stat package
    stats_pp=stats::PP.test(rnorm(100))
    
    
    
    
    # extract p value from different test result. 
    p_value_vector=c(test_fractal=attr(test_fractal,"pvals")[1],tseries_adf=tseries_adf$p.value,
                     tseries_kpss=tseries_kpss$p.value,tseries_pp=tseries_pp$p.value,
                     stats_pp=stats_pp$p.value)
    return(p_value_vector)
  }
  
  colnames(p_value_matrix)=c("spectra_method","tseries_adf","tseries_kpss","tseries_pp","stats_pp")
  
  return(p_value_matrix)
}

P_value_matrix()
aTSA::adf.test(rnorm(100))
