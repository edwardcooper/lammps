## write the parallel version of map function map_c with foreach and doParallel. map_c means along column. 
map_pc=function(data,FUN,nthread=detectCores(),type="FORK"){
  library(foreach)
  library(doParallel)
  
  # make a fork cluster with half of the cores.
  # Half of the cores are the number of real cores. 
  cl=makeCluster(nthread,type=type)
  registerDoParallel(cl)
  
  #result should return a list by default.
  out_list=foreach(i=1:ncol(data))%dopar%{
    result=FUN(data[,i])
  }
  message=paste("The computer is using ", getDoParWorkers(), " core(s).",sep="")
  print(message)
  # close the parallel ports 
  gc()
  stopCluster(cl)
  stopImplicitCluster()
  
  return(out_list)
}
# example use 
# random_data=data.frame(matrix(rnorm(1e8),ncol=10,nrow=1e7))
#                                     
# cumsum_result=map_p(data=random_data,sum)
# 
# cumsum_result%>%str



