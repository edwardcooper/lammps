## Use this script to calculate the power of different stationary test. 


# I will vary the size of the sample and number of tests to see how they affect the power of statistical test. 

# due to the limit on the p value produced the type I error we could choose from is (0.01,0.1)

# write the p_value of statistical test into a function P_value_matrix(). 

P_value_rnorm=function(num_test=10,num_sample=1000){
  library(foreach)
 
  # generate the random number and do the test in the same loop. Reduce memory usage. 
  p_value_matrix=foreach(i=1:num_test,.combine=rbind,.packages = c("fractal","tseries"))%do%{
    
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
    stats_pp=stats::PP.test(random_number)
    
    
    
    # extract p value from different test result. 
    p_value_vector=c(test_fractal=attr(test_fractal,"pvals")[1],tseries_adf=tseries_adf$p.value,
                     tseries_kpss=tseries_kpss$p.value,tseries_pp=tseries_pp$p.value,
                     stats_pp=stats_pp$p.value)
    
    return(p_value_vector)
  }
  
  #####################################################################################
  #information about the real process, null and alternative hypothesis for each test.
  true_value=c("null","alternative","null","alternative","alternative")
  null=c("stationary","unit root present","level stationary","unit root present","unit root present")
  alternative=c("non-stationary","stationary","unit root present","stationary","stationary")
  sample_size=rep(num_sample,5) # five tests 
  num_test=rep(num_test,5) # five tests
  #####################################################################################
  
  p_value_matrix=rbind(true_value,null,alternative,sample_size,num_test,p_value_matrix)
  colnames(p_value_matrix)=c("spectra_method","tseries_adf","tseries_kpss","tseries_pp","stats_pp")
  
  return(p_value_matrix)
}

#data=P_value_rnorm()

# I leave out the methods from aTSA package for future reference. 
# I leave out urca package since it took a very long time to compute and only outputs a test statistic. 
# I leave out locits package since this method has some length restrictions.
## The length of time series must be a number that is 2^n for locits package. The author say he will do it for arbitrary length in the future. 




# Extend the scope of the test. Different distributions and ARMA models. 

















# This function reads the p value matrix, null hypothesis, and alternative hypothesis, 
# then it returns the power of the statistical test, 
# and the result reject the null hypothesis or fail to reject the null hypothesis.



# alpha is set between 0.01 and 0.1. Restricted by the outcome of statistical test.
# we will calculate the power of test from the true_value and p values generated. 

power_of_test=function(data,alpha=0.05){
  
  # extract information about statistical test from the data
  information=data[1:5,]
  
  
  # calculate power of test 
  ## PSR test. A.K.A spectra method
  spectra_method=data[,"spectra_method"]>alpha
  ## adf test
  tseries_adf=data[,"tseries_adf"]<alpha
  ## kpss test 
  tseries_kpss=data[,"tseries_kpss"]>alpha
  ## pp test 
  tseries_pp=data[,"tseries_pp"]<alpha
  ## stats pp test
  stats_pp=data[,"stats_pp"]<alpha
  
  # combine the above result into a data.frame. 
  test_logic_data=data.frame(spectra_method=spectra_method,tseries_adf=tseries_adf,
                             tseries_kpss=tseries_kpss,tseries_pp=tseries_pp,stats_pp=stats_pp)
  
  # write a function to calculate the percentage of tests that does it correctly. 
  true_percent=function(data){
    return( sum(data)/length(data) )
  }
  
  power_of_test=apply(test_logic_data,2,true_percent)
  result=rbind(information, power_of_test)
  colnames(result)=colnames(data)
  return(result)
}

# Extend the functionality of power_of_test to use ture_value and p value to get the power of test. 
# The function here is not exactly a power of test, which is 1-beta, but a test of different tests to distinguish the stationarity. 
# library(magrittr)
# 
# power_test=P_value_rnorm(num_test=1000,num_sample=1000)%>%power_of_test(alpha=0.02)
# 
# power_test
