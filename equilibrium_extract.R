# The purpose of this function is to get the part of simulation that is in equilibrium from the energies, temperatures and etc. 

# Could only chosee alpha in (0.01, 0.1)
equilibrium_extract=function(data, del_num=10,alpha=0.05){
  # do the statistical test on the data
  test_fractal=fractal::stationarity(data) #null hypothesis: stationary. 
  # extract the p value. 
  p_value=attr(test_fractal,"pvals")[1]
  # test to see if the p value is smaller than alpha (type I error rate). 
  # If it is TRUE, it means that the p value is smaller than alpha, then null hypothesis that the data is staionary is rejected. Data not stationary. 
  # If it is FALSE, it means that the p value is bigger than alpha, then we fail to reject the null hypothesis that the data is staionary. Data stationary. 
  
  
  non_stationary=p_value<alpha
  
  # If it is stationary, then we will return the original data.
  # If it is not stationary, then we will need to chop off some of the data to test again until we find the part of data that is stationary. 
  if(!non_stationary){
    
    result_data=data
    print("This data is stationary")
  }else if (non_stationary){
    
    # while non_stationary is 
    while(non_stationary){
      del_index=1:del_num
      data=data[-del_num]
      test_fractal=fractal::stationarity(data)
      p_value=attr(test_fractal,"pvals")[1]
      non_stationary=p_value>alpha
    }
    
    result_data=data
    print("The length of the data is", length(data))
  }
  gc()
  return(result_data)
}


rnorm(100)%>%equilibrium_extract()

