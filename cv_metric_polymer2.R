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
 cv_metric_polymer2=function(
  Path="~/Dropbox/lammps"
  ,polymer="PMMA_big"
  ,temp=seq(300,600,by=20)
  ,filename="MSD.g0.colmean.Tmean.1.txt"
   ,method="bs_splines"
  ,timesteps=500:4001
  ,k=10
 , target="variance_ratio_linear" 
  ){
  MSD.PS.g0=Data.import(Path=Path, filename = filename, temp=temp,polymer=polymer)
  # subset the dataset to contain only data for those timestep
  MSD.PS.g0=MSD.PS.g0[timesteps,]%>%select(-time)
  ################################################################
  
  # 
  MSD.PS.g0=cbind(MSD.PS.g0,time_steps=timesteps)
  # change the time units to ps
  # change  the length unit to nm.
  library(foreach)
  colnames(MSD.PS.g0)=c(paste("T",temp,sep=""),"time_steps")
  
 
  
  
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
                              ,rousse_einstein={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])]   ) }
                              ,rousse_reptation={ lm(train_data[,temperature]~I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,reptation_einstein={ lm(train_data[,temperature]~I(time_steps^(1/4))+I(time_steps),data=train_data[,c(temperature,dim(train_data)[2])] ) }
                              ,hinges={library(earth);earth(train_data[,temperature]~I(time_steps)+I(time_steps^(1/2))+I(time_steps^(1/4)),data=train_data[,c(temperature,dim(train_data)[2])] )}
                              ,ns_splines={library(splines);lm(train_data[,temperature]~ns(I(time_steps^(1/2))+I(time_steps^(1/4))+I(time_steps)),data=train_data[,c(temperature,dim(train_data)[2])])}
                              ,bs_splines={library(splines);lm(train_data[,temperature]~bs(I(time_steps^(1/2))+I(time_steps^(1/4))+I(time_steps)),data=train_data[,c(temperature,dim(train_data)[2])])}
      )
      
      library(magrittr)
      paste("temperature:",temperature,"Fold:",fold,sep="")%>%message
      
      
      
     
 
      r_squared=1-mean((predict(regression_model,train_data)-train_data[,temperature])^2)/mean((mean(train_data[,temperature])-train_data[,temperature])^2)
      prediction_msd=predict(regression_model,newdata=test_data)
      
      
      predicted_r_squared=1-mean((prediction_msd-test_data[,temperature])^2)/mean(( test_data[,temperature]-mean(test_data[,temperature]) )^2)
      
      # all computation to do squared error comparison. 
      
      # linear models performance. aka Einstein model 
      regression_model_linear=lm(train_data[,temperature]~I(time_steps),data=train_data )
      prediction_msd_linear=predict(regression_model_linear,newdata=test_data%>%select(time_steps))
      variance_ratio_einstein=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_linear-test_data[,temperature])^2)
      # Rousse model performance
      regression_model_rousse=lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data )
      prediction_msd_rousse=predict(regression_model_rousse,newdata=test_data%>%select(time_steps))
      variance_ratio_rousse=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_rousse-test_data[,temperature])^2)
      # reptation model performance 
      regression_model_reptation=lm(train_data[,temperature]~I(time_steps^(1/2)),data=train_data )
      prediction_msd_reptation=predict(regression_model_reptation,newdata=test_data%>%select(time_steps))
      variance_ratio_reptation=sum((prediction_msd-test_data[,temperature])^2) /sum((prediction_msd_reptation-test_data[,temperature])^2)
      
      cv_metrics=c(r_squared,predicted_r_squared
                   ,variance_ratio_einstein,variance_ratio_rousse,variance_ratio_reptation
                   ,temp[temperature])
   }
 
  }

  
  
  rownames(all_temperatures_cv_metrics)=1:dim(all_temperatures_cv_metrics)[1]
  
  
  # cross-validation loop
  
  
  
  all_temperatures_cv_metrics=all_temperatures_cv_metrics%>%as.data.frame()
  colnames(all_temperatures_cv_metrics)=c("r_squared","predicted_r_squared"
                                          ,"variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation"
                                          ,"temperature")
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

# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="bs_splines",timesteps=500:4001,k=10,target=c("r_squared","predicted_r_squared"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="bs_splines",timesteps=500:4001,k=10,target=c("variance_ratio_reptation","variance_ratio_linear"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                             ,method="ns_splines",timesteps=500:4001,k=10,target=c("r_squared","predicted_r_squared"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="ns_splines",timesteps=500:4001,k=10,target=c("variance_ratio_reptation","variance_ratio_linear"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="hinges",timesteps=500:4001,k=10,target=c("r_squared","predicted_r_squared"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="hinges",timesteps=500:4001,k=10,target=c("variance_ratio_reptation","variance_ratio_linear"))
# 
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="mixall",timesteps=500:4001,k=10,target=c("r_squared","predicted_r_squared"))
# 
# cv_metric_polymer2(Path="~/Dropbox/lammps",polymer="PMMA_big", temp=seq(300,600,by=20),filename="MSD.g0.colmean.Tmean.1.txt"
#                    ,method="mixall",timesteps=500:4001,k=10,target=c("variance_ratio_reptation","variance_ratio_linear"))

four_combined_plots=function(method="mixall"
                            ,timesteps=200:4001
                            ,polymer="PS_20"
                            ,temp=seq(200,500,by=20)
                            ,k=10  
                            ,variance_variables=c("variance_ratio_einstein","variance_ratio_rousse","variance_ratio_reptation")
){
  
  
  pmma_r_squared_g0=cv_metric_polymer2(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g0.colmean.Tmean.1.txt"
                                      ,method=method,timesteps=timesteps,k=k,target=c("r_squared","predicted_r_squared") )
  
  pmma_variance_ratio_g0=cv_metric_polymer2(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g0.colmean.Tmean.1.txt"
                                           ,method=method,timesteps=timesteps,k=k,target=variance_variables )
  

  
  pmma_r_squared_g1=cv_metric_polymer2(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g1.colmean.Tmean.1.txt"
                                      ,method=method,timesteps=timesteps,k=k,target=c("r_squared","predicted_r_squared") )
  
  pmma_variance_ratio_g1=cv_metric_polymer2(Path="~/Dropbox/lammps",polymer=polymer, temp=temp,filename="MSD.g1.colmean.Tmean.1.txt"
                                           ,method=method,timesteps=timesteps,k=k,target=variance_variables )
  
  
  
  
  
  
  
  # library(latex2exp)
  r_squared_pmma_g0_all=rbind(ggplotGrob(pmma_r_squared_g0), ggplotGrob(pmma_variance_ratio_g0)
                              , size = "first")
  
  r_squared_pmma_g1_all=rbind(ggplotGrob(pmma_r_squared_g1), ggplotGrob(pmma_variance_ratio_g1)
                              , size = "first")
  
  library(gridExtra)
  final_combined_plot=grid.arrange(r_squared_pmma_g0_all,r_squared_pmma_g1_all,nrow=1)
  
  
  return(final_combined_plot)
  
  
}


# 
# four_combined_plots(method="mixall",timesteps=200:4001  ,polymer="PS_20",temp=seq(200,500,by=20) ,k=5  )
#                   
# 
# 
# four_combined_plots(method="hinges",timesteps=200:4001  ,polymer="PS_20",temp=seq(200,500,by=20) ,k=5 )
#                                   
# 
# four_combined_plots(method="ns_splines",timesteps=200:4001 ,polymer="PS_20",temp=seq(200,500,by=20) ,k=5 )
# 
# 
# four_combined_plots(method="bs_splines",timesteps=200:4001  ,polymer="PS_20",temp=seq(200,500,by=20) ,k=5 )
                   