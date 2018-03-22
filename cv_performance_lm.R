# define a function to read in data
library(foreach)
library(dplyr)
Data.import=function(Path="~/Dropbox/lammps" , filename = "MSD.colmean.matrix.1.txt", temp=c(320,345,370,395,420,445,470,500,550,600),polymer="PMMA"){
  
  # read in data into a matrix
  MSDcol=foreach(i=1:length(temp),.combine = cbind)%do% {
    # set correct path for the data file
    path=paste(Path,"/", polymer,"/atom",temp[i], sep='')
    path
    setwd(path)
    # read in the data
    MSDcol=read.table(file=filename,header=TRUE,sep=',')
    library(magrittr)
    MSDcol%>%dim
    MSDcol=MSDcol[[1]]
    return(MSDcol)
  }
  
  # add a time variable 
  time=1:dim(MSDcol)[1]
  # change data into data frame 
  MSDcol=cbind(time,MSDcol)%>%as.data.frame()
  
}



cv_r_squared=function(filename1="MSD.g0.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="mixdiff",
                      timesteps=1200:5001,k=20){
  MSD.PS.g0=Data.import(filename=filename1, temp=temp,polymer=polymer)
  
  MSD.PS.g0=MSD.PS.g0%>%as.data.frame()
  
  # subset the dataset to contain only data for those timestep
  MSD.PS.g0=MSD.PS.g0[timesteps,]%>%select(-time)
  ################################################################
  
  
  MSD.PS.g0=cbind(MSD.PS.g0,time_steps=seq_along(timesteps))
  # change the time units to ps
  # change  the length unit to nm. 
  MSD.PS.g0=MSD.PS.g0/100
  library(foreach)
  
  index_list=caret::createFolds(timesteps, k = k, list = TRUE, returnTrain = FALSE)
  
  # a for loop here for cross-validation. Use for loop for this or foreach? 
  # a for loop here for different temperature. This is the outer loop. 
  
  
  # different temperature loop
  
  all_temperatures_cv_metrics=foreach(temperature=seq_along(temp)+1,.combine = rbind)%do%{
    
    
    cv_performance_metrics=foreach(fold=1:k,.combine = rbind)%do%{
      
      index_fold=index_list[[fold]]
      test_data=MSD.PS.g0[index_fold,]
      train_data=MSD.PS.g0[-index_fold,]
      
      
      regression_model=switch(method,
                              linear={ lm(train_data[,temperature]~-1+I(time_steps^2)+I(time_steps^4)+I(time_steps^6)+I(time_steps^8)+I(time_steps^10),data=train_data ) }
                              ,frac={ lm(train_data[,temperature]~-1+I(time_steps^(3/2))+I(time_steps^(7/2))+I(time_steps^(11/2))+I(time_steps^(15/2))+I(time_steps^(19/2)),data=train_data ) }
                              ,subdiff={ lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data ) }
                              ,diff={ lm(train_data[,temperature]~I(time_steps),data=train_data ) }
                              ,mixdiff={ lm(train_data[,temperature]~I(time_steps)+I(time_steps^(1/2)),data=train_data ) }
      )
      
      
      
      
      r_squared=regression_model%>%summary%>%.$r.squared
      adj_r_squared=regression_model%>%summary%>%.$adj.r.squared
      
      prediction_msd=predict(regression_model,newdata=test_data%>%select(time_steps))
      
      predicted_r_squared=1-sum((prediction_msd-test_data[,temperature])^2)/sum((test_data[,temperature]-mean(test_data[,temperature]))^2)
      
      cv_metrics=c(r_squared,adj_r_squared,predicted_r_squared,temp[temperature-1])
      
      return(cv_metrics)
    }
    
  }
  
  colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared","temperature")
  rownames(all_temperatures_cv_metrics)=1:dim(all_temperatures_cv_metrics)[1]
  # cross-validation loop 
  library(ggplot2)
  # all_temperatures_cv_metrics=all_temperatures_cv_metrics%>%as.data.frame()
  
  all_temperatures_cv_metrics=all_temperatures_cv_metrics%>%as.data.frame()%>%select(temperature,predicted_r_squared)
  all_temperatures_cv_metrics$temperature=as.factor(all_temperatures_cv_metrics$temperature)
  
 # all_temperatures_cv_metrics%>%ggplot()+geom_boxplot(aes(x=temperature,y=predicted_r_squared,group=temperature),color="red")
  
 
  # print(all_temperatures_cv_metrics)
  library(ggplot2)
  
  # graph=ggplot(data=all_temperatures_cv_metrics%>%reshape2::melt(id.vars="temperature"))+geom_boxplot(aes(x=temperature,y=value,fill=variable))+
  #   scale_y_continuous(limits = c(-1,1))
  graph=ggplot(data=all_temperatures_cv_metrics)+geom_boxplot(aes(x=temperature,y=predicted_r_squared),fill="red")+xlab("Temperature(K)")
  
  
  return(graph)
  
}



# comparing different models. 







## PMMA 

#360K
r_squared_pmma_g0_subdiff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="subdiff",timesteps=1000:5001,k=10)
#360K
r_squared_pmma_g0_diff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="diff",timesteps=1000:5001,k=10)
#360K
r_squared_pmma_g0_mixdiff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="mixdiff",timesteps=1000:5001,k=10)

# # No clear conclusion
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="subdiff",timesteps=1000:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="diff",timesteps=1200:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="mixdiff",timesteps=1200:5001,k=10)#+ scale_y_continuous(limits = c(0.9,1))
# 


# # 400K
# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="subdiff",timesteps=1200:5001,k=10)
# #420K
# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="diff",timesteps=1200:5001,k=10)
# #380K
# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="mixdiff",timesteps=1200:5001,k=10)


# # 360K 
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="subdiff",timesteps=1200:5001,k=10)# + scale_y_continuous(limits = c(0.75,1))
# # 360K 
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="diff",timesteps=1200:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))
# # 360K 
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(300,600,by=20),polymer="PMMA_big",method="mixdiff",timesteps=1200:5001,k=10)# + scale_y_continuous(limits = c(0.75,1))
# 

## PS


# 320K +440K
r_squared_ps_g0_subdiff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="subdiff",timesteps=1000:5001,k=10)

# 340K +440K
r_squared_ps_g0_diff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="diff",timesteps=1000:5001,k=10)
# 320K +440K
r_squared_ps_g0_mixdiff=cv_r_squared(filename1="MSD.g0.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="mixdiff",timesteps=1000:5001,k=10)

# # 280K,320K
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="subdiff",timesteps=1000:5001,k=10)#+ scale_y_continuous(limits = c(0,1))
# # 280K, 320K
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="diff",timesteps=1000:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))
# # 280k, 320K, 380K
# cv_r_squared(filename1="MSD.g0.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="mixdiff",timesteps=1200:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))



# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="subdiff",timesteps=1000:5001,k=10)
# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="diff",timesteps=1000:5001,k=10)
# cv_r_squared(filename1="MSD.g1.colmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="mixdiff",timesteps=1000:5001,k=10)

# # 320K
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="subdiff",timesteps=1000:5001,k=10)# + scale_y_continuous(limits = c(0.75,1))
# # 320K
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="diff",timesteps=1000:5001,k=10)#+ scale_y_continuous(limits = c(0.75,1))
# # 320K, 380K
# cv_r_squared(filename1="MSD.g1.colmean.Tmean.1.txt",temp=seq(200,500,by=20),polymer="PS_20",method="mixdiff",timesteps=1000:5001,k=10)# + scale_y_continuous(limits = c(0.75,1))



r_squared_pmma_all=rbind(ggplotGrob(r_squared_pmma_g0_subdiff), ggplotGrob(r_squared_pmma_g0_diff),ggplotGrob(r_squared_pmma_g0_mixdiff), size = "first")

r_squared_ps_all=rbind(ggplotGrob(r_squared_ps_g0_subdiff), ggplotGrob(r_squared_ps_g0_diff),ggplotGrob(r_squared_ps_g0_mixdiff), size = "first")




# PS 
grid.arrange(r_squared_pmma_all,r_squared_ps_all,nrow=1)



