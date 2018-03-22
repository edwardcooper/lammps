

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

#====================================================================================================


cv_metric_polymer=function(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
                           ,method="mixall",timesteps=500:4001,k=10,target ){
  MSD.PS.g0=Data.import(Path=Path, filename = filename, temp=temp,polymer=polymer)
  # subset the dataset to contain only data for those timestep
  MSD.PS.g0=MSD.PS.g0[timesteps,]%>%select(-time)
  ################################################################
  
  # 
  MSD.PS.g0=cbind(MSD.PS.g0,time_steps=timesteps)
  # change the time units to ps
  # change  the length unit to nm.
  MSD.PS.g0=MSD.PS.g0
  library(foreach)
  colnames(MSD.PS.g0)=c(paste("T",temp,sep=""),"time_steps")
  
  # #regular plot 
  # ###########################################
  # 
  # MSD.PS.g0.melt=reshape2::melt(MSD.PS.g0,id.vars="time_steps")
  # library(dplyr)
  # 
  # library(ggplot2)
  # ggplot(data=MSD.PS.g0.melt,aes(x=time_steps,y=value,colour=variable)) +geom_point(size=0.5)+
  #   ylab("MSD")+
  #   ggtitle("PS g0 MSD vs time plot")
  
  
  index_list=caret::createFolds(timesteps, k = k, list = TRUE, returnTrain = FALSE)
  
  # this function assumes that the time_steps are at the end of the dataset. 
  
  all_temperatures_cv_metrics=foreach(temperature=seq_along(temp),.combine = rbind)%do%{
    
    
    cv_performance_metrics=foreach(fold=1:k,.combine = rbind)%do%{
      
      index_fold=index_list[[fold]]
      test_data=MSD.PS.g0[index_fold,]
      train_data=MSD.PS.g0[-index_fold,]
      
      # the order of the term is reflected in the name. 
      regression_model=switch(method,
                              mixall={ lm(train_data[,temperature]~I(time_steps)+I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse={ lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,reptation={ lm(train_data[,temperature]~I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,einstein={ lm(train_data[,temperature]~I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse_einstein={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse_reptation={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,reptation_einstein={ lm(train_data[,temperature]~I(time_steps^(1/4))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
      )
      
      library(magrittr)
      paste("temperature:",temperature,"Fold:",fold,sep="")%>%message
      
      r_squared=regression_model%>%summary%>%.$r.squared
      adj_r_squared=regression_model%>%summary%>%.$adj.r.squared
      
      
      coefficients_from_model=switch(method
                                     ,mixall={ c(regression_model$coefficients[2],regression_model$coefficients[3],regression_model$coefficients[4]) }
                                     ,rousse={ c(regression_model$coefficients[2]) }
                                     ,reptation={ c(regression_model$coefficients[2])  }
                                     ,einstein={ c(regression_model$coefficients[2]) }
                                     ,rousse_einstein={ c(regression_model$coefficients[2],regression_model$coefficients[3]) }
                                     ,rousse_reptation={ c(regression_model$coefficients[2],regression_model$coefficients[3]) }
                                     ,reptation_einstein={ c(regression_model$coefficients[2],regression_model$coefficients[3]) }
      )
      
      
      prediction_msd=predict(regression_model,newdata=test_data%>%select(time_steps))
      predicted_r_squared=1-mean((prediction_msd-test_data[,temperature])^2)/mean(( test_data[,temperature]-mean(test_data[,temperature]) )^2)
      
      
      # all computation to do squared error comparison. 
      
      # linear models performance. aka Einstein model 
      regression_model_linear=lm(train_data[,temperature]~I(time_steps),data=train_data )
      prediction_msd_linear=predict(regression_model_linear,newdata=test_data%>%select(time_steps))
      variance_ratio_linear=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_linear-test_data[,temperature])^2)
      # Rousse model performance
      regression_model_rousse=lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data )
      prediction_msd_rousse=predict(regression_model_rousse,newdata=test_data%>%select(time_steps))
      variance_ratio_rousse=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_rousse-test_data[,temperature])^2)
      # reptation model performance 
      regression_model_reptation=lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data )
      prediction_msd_reptation=predict(regression_model_reptation,newdata=test_data%>%select(time_steps))
      variance_ratio_reptation=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_reptation-test_data[,temperature])^2)
      
      cv_metrics=c(r_squared,adj_r_squared,predicted_r_squared
                   ,variance_ratio_linear,variance_ratio_rousse,variance_ratio_reptation
                   ,coefficients_from_model,temp[temperature])
      
      return(cv_metrics)
    }
    
  }
  
  # 
  if(method=="mixall"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                ,"coefficients_from_model_einstein","coefficients_from_model_rousse","coefficients_from_model_reptation"
                                                                ,"temperature")}
  if(method=="rousse"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                ,"coefficients_from_model_rousse"
                                                                ,"temperature")}
  if(method=="reptation"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                   ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                   ,"coefficients_from_model_reptation"
                                                                   ,"temperature")}
  if(method=="einstein"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                  ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                  ,"coefficients_from_model_einstein"
                                                                  ,"temperature")}
  if(method=="rousse_einstein"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                         ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                         ,"coefficients_from_model_rousse","coefficients_from_model_einstein"
                                                                         ,"temperature")}
  if(method=="rousse_reptation"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                          ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                          ,"coefficients_from_model_rousse","coefficients_from_model_reptation"
                                                                          ,"temperature")}
  if(method=="reptation_einstein"){ colnames(all_temperatures_cv_metrics)=c("r_squared","adj_r_squared","predicted_r_squared"
                                                                            ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                                                            ,"coefficients_from_model_reptation","coefficients_from_model_einstein"
                                                                            ,"temperature")}
  
  rownames(all_temperatures_cv_metrics)=1:dim(all_temperatures_cv_metrics)[1]
  
  
  # cross-validation loop
  

  
  all_temperatures_cv_metrics=all_temperatures_cv_metrics%>%as.data.frame()
  all_temperatures_cv_metrics=all_temperatures_cv_metrics[,c("temperature",target)]
  
  
  all_temperatures_cv_metrics.melt=reshape2::melt(all_temperatures_cv_metrics,id.vars="temperature")
  
  all_temperatures_cv_metrics.melt$temperature=as.factor(all_temperatures_cv_metrics.melt$temperature)
  
  # return(all_temperatures_cv_metrics.melt)
  library(ggplot2)
  
  
  graph=ggplot(data=all_temperatures_cv_metrics.melt,aes(x=temperature,y=value,colour=variable))+geom_boxplot()+
                                                                          xlab("Temperature(K)")+
                                                                          theme(text = element_text(size=10))
  
  return(graph)
}




six_combined_plots=function(method="mixall"
                            ,timesteps=200:4001
                            ,polymer="PS_20"
                            ,temp=seq(200,500,by=20)
                            ,k=10  
                            ,coefficient_variables=c("coefficients_from_model_einstein","coefficients_from_model_rousse","coefficients_from_model_reptation")
                            ,variance_variables=c("variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation")
                            ){
  
  
  pmma_r_squared_g0=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g0.colmean.Tmean.1.txt"
                                      ,method=method,timesteps=timesteps,k=k,target=c("r_squared","adj_r_squared","predicted_r_squared") )
  
  pmma_variance_ratio_g0=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g0.colmean.Tmean.1.txt"
                                           ,method=method,timesteps=timesteps,k=k,target=variance_variables )
  
  pmma_coefficients_g0=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g0.colmean.Tmean.1.txt"
                                         ,method=method,timesteps=timesteps,k=k,target=coefficient_variables )
  
  pmma_r_squared_g1=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g1.colmean.Tmean.1.txt"
                                      ,method=method,timesteps=timesteps,k=k,target=c("r_squared","adj_r_squared","predicted_r_squared") )
  
  pmma_variance_ratio_g1=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g1.colmean.Tmean.1.txt"
                                           ,method=method,timesteps=timesteps,k=k,target=variance_variables )
  
  pmma_coefficients_g1=cv_metric_polymer(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g1.colmean.Tmean.1.txt"
                                         ,method=method,timesteps=timesteps,k=k,target=coefficient_variables )
  
  
  
  
  
  
  # library(latex2exp)
  r_squared_pmma_g0_all=rbind(ggplotGrob(pmma_r_squared_g0), ggplotGrob(pmma_variance_ratio_g0)
                              ,ggplotGrob(pmma_coefficients_g0)
                              , size = "first")
  
  r_squared_pmma_g1_all=rbind(ggplotGrob(pmma_r_squared_g1), ggplotGrob(pmma_variance_ratio_g1)
                              ,ggplotGrob(pmma_coefficients_g1)
                              , size = "first")
  
  library(gridExtra)
  final_combined_plot=grid.arrange(r_squared_pmma_g0_all,r_squared_pmma_g1_all,nrow=1)
  
  
  return(final_combined_plot)
  
  
}


# methods=c("mixall","rousse","reptation","einstein","rousse_einstein","rousse_reptation","reptation_einstein")
# change the coefficient_variables for different models you want to use. 
# example use 
# coefficient_variables=c("coefficients_from_model_einstein"
#                         ,"coefficients_from_model_rousse"
#                         ,"coefficients_from_model_reptation")
# 
# variance_variables=c("variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation")
# 
# mix_all_ps=six_combined_plots(method="mixall"
#                               ,timesteps=200:4001
#                               ,polymer="PS_20"
#                               ,temp=seq(200,500,by=20)
#                               ,k=10  
#                               ,coefficient_variables=coefficient_variables
#                               ,variance_variables=variance_variables)
# mix_all_ps

# cv predicted values plotting function. 

cv_predict_plot=function(data,newdata){
  data=data%>%as.data.frame()
  data=data%>%filter(time>1)
  # add correct colname 
  temp=seq(300,600,by=20)
  colnames(data)=c("time",paste("T",temp,sep=""))
  
  
  #regular plot 
  ###########################################
  library(reshape2)
  data.melt=melt(data,id.vars="time")
  newdata.melt=melt(newdata,id.vars = "time_steps")
  
  library(ggplot2)
  graph=ggplot(data=data.melt,aes(x=(time),y=(value),colour=variable)) +geom_point(size=0.05)+
    geom_line(data=newdata.melt,aes(x=(time_steps),y=(value),color=variable))+
    ylab("MSD")+
    ggtitle("PMMA g0 Tmean MSD vs time plot")
  return(graph)
}



cv_predict_polymer=function(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
                            ,method="mixall",timesteps=500:4001,k=10 ){
  MSD.PS.g0=Data.import(Path=Path, filename = filename, temp=temp,polymer=polymer)
  # subset the dataset to contain only data for those timestep
  MSD.PS.g0=MSD.PS.g0[timesteps,]%>%select(-time)
  ################################################################
  
  # 
  MSD.PS.g0=cbind(MSD.PS.g0,time_steps=timesteps)
  # change the time units to ps
  # change  the length unit to nm.
  MSD.PS.g0=MSD.PS.g0
  library(foreach)
  colnames(MSD.PS.g0)=c(paste("T",temp,sep=""),"time_steps")
  
  # #regular plot 
  # ###########################################
  # 
  # MSD.PS.g0.melt=reshape2::melt(MSD.PS.g0,id.vars="time_steps")
  # library(dplyr)
  # 
  # library(ggplot2)
  # ggplot(data=MSD.PS.g0.melt,aes(x=time_steps,y=value,colour=variable)) +geom_line(size=0.5)+
  #   ylab("MSD")+
  #   ggtitle("PS g0 MSD vs time plot")
  
  
  index_list=caret::createFolds(timesteps, k = k, list = TRUE, returnTrain = FALSE)
  
  # this function assumes that the time_steps are at the end of the dataset. 
  
  all_temperatures_cv_metrics=foreach(temperature=seq_along(temp),.combine = cbind)%do%{
    
    
    cv_performance_metrics=foreach(fold=1:k,.combine = rbind)%do%{
      
      index_fold=index_list[[fold]]
      test_data=MSD.PS.g0[index_fold,]
      train_data=MSD.PS.g0[-index_fold,]
      
      # the order of the term is reflected in the name. 
      regression_model=switch(method,
                              mixall={ lm(train_data[,temperature]~I(time_steps)+I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse={ lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,reptation={ lm(train_data[,temperature]~I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,einstein={ lm(train_data[,temperature]~I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse_einstein={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,rousse_reptation={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,reptation_einstein={ lm(train_data[,temperature]~I(time_steps^(1/4))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
      )
      
      library(magrittr)
      paste("temperature:",temperature,"Fold:",fold,sep="")%>%message
      
      
      paste("The number of data in test data:", dim(test_data)[1])%>%message
      
      prediction_msd=predict(regression_model,newdata=test_data%>%select(time_steps))
      predicted_r_squared=1-mean((prediction_msd-test_data[,temperature])^2)/mean(( test_data[,temperature]-mean(test_data[,temperature]) )^2)
      
      #msd_name=paste("T",temp[temperature],sep="")
      
      prediction_msd_with_timesteps=cbind(time_steps=test_data%>%select(time_steps),msd=prediction_msd)
      
      return(prediction_msd_with_timesteps)
    }
    # reorder the predicted value by time_steps
    cv_performance_metrics=cv_performance_metrics%>%arrange(time_steps)
    # get rid of the time_steps variable for easier combination of data. 
    cv_performance_metrics=cv_performance_metrics%>%select(-time_steps)
  }
  
  # add time_steps variable after combining all temperature results. 
  all_temperatures_cv_metrics=cbind(all_temperatures_cv_metrics,time_steps=timesteps)
  colnames(all_temperatures_cv_metrics)=c(paste("T",temp,sep=""),"time_steps")
  
  # return value of the entire function.   
  return(all_temperatures_cv_metrics)
}

