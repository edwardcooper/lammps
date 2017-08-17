# read in data from the log file to examine the results 

tvp_data=function(file="tvp.csv"){
  tvp=read.table(file)
  
  
  colnames(tvp)=c("property","value")
  
  library(dplyr)
  t=tvp%>%filter(property=="Temp")
  
  
  
  p=tvp%>%filter(property=="Press")
  
  v=tvp%>%filter(property=="Volume")
  
  table(tvp$property)%>%print()
  
  tvp=cbind(timestep=1:(dim(t)[1]), temperature=t$value, pressure=p$value,volume=v$value)
  
  tvp=tvp%>%as.data.frame()
  
  return(tvp)
  }

tvp_plot=function(file="tvp.csv"){
  tvp=tvp_data(file)
  
  library(ggplot2)
  
  p1=tvp%>%ggplot()+geom_line(aes(x=timestep,y=temperature),color="blue")
  p2=tvp%>%ggplot()+geom_line(aes(x=timestep,y=pressure),color="red")
  
  
  library(grid)
  grid.newpage()
  plot=grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"))
  return(plot)
}

tvp_plot(file="tvp320solomon.csv")
tvp_plot(file="tvp345solomon.csv")
tvp_plot(file="tvp320Sbond.csv")
tvp_plot(file="tvp320Rbond.csv")

# data=tvp_data(file="tvp320Rbond.csv")
# head(data)
