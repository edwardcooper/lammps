# Testing som codes to see the performance for optimization. 


setwd("~/Dropbox/lammps/PMMA_big/atom600")
source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")

MSD=function(data){
  (data$xu-data$xu[1])^2+ (data$yu-data$yu[1])^2+ (data$zu-data$zu[1])^2
}
library(foreach)
library(magrittr)
library(microbenchmark)






timeRecordB()
### Load the data
atom.600.1=read.table(file="atom.600_1.txt", header=FALSE, sep=" ", stringsAsFactors = FALSE, fill=TRUE
                      ,col.names = c("atom-id","type","mol","xu","yu","zu")
                      ,colClasses = c("numeric","numeric","numeric","numeric","numeric","numeric")
                      ,na.strings = c("ITEM:","TIMESTEP","NUMBER","OF","ATOMS","BOX","BOUNDS", "pp","id","type","mol","xu","yu","zu")
)
timeRecordB(output_message = "Load in data")
# Clear out all rows with NAs.
atom.600.1=atom.600.1[complete.cases(atom.600.1),]
summary(atom.600.1)
gc()
timeRecordR(unit="min")%>%select(output_message,run_time)

atom.600.1%>%sapply(class)
library(microbenchmark)

microbenchmark(atom.600.1%>%filter(atom.id==1)%>%select(xu,yu,zu)%>%MSD
               ,atom.600.1[atom.600.1[,"atom.id"]==1,]%>%.[,c("xu","yu","zu")]%>%MSD
               ,atom.600.1[atom.600.1[,"atom.id"]==1,]%>%select(xu,yu,zu)%>%MSD
               )
# use the first expression, which is the fastest. 


# ## This is faster than select
# atom.600.1[,"atom.id"]%>%max()
# atom.600.1[,"atom.id"]%>%min()
# atom.600.1[,"mol"]%>%max()
# atom.600.1[,"mol"]%>%min()



microbenchmark(
  
  MSD.all.matrix.600=foreach(n=1:2,.combine=rbind)%do%{ 
    atom.600.1%>%filter(atom.id==n)%>%select(xu,yu,zu)%>%MSD
  } 
  ,
  MSD.all.matrix.600.2=foreach(n=1:2,.combine=rbind)%do%{ 
    atom.600.1%>%filter(atom.id==n)%>%select(xu,yu,zu)
  }
  ,times = 10
)
# These two expressions use about the same amount of time. Thus, 

## The performance bottleneck is not calculating the MSD but searching for all atom.id==1 and select xu,yu,zu values to calculate MSD 


## Try to use data.table instead of dplyr 
library(data.table)
atom.600.1.data_table=atom.600.1%>%as.data.table
timeRecordB()
microbenchmark(
  
  MSD.all.matrix.600=foreach(n=1:10,.combine=rbind)%do%{ 
    atom.600.1%>%filter(atom.id==n)%>%select(xu,yu,zu)%>%MSD
  } 
  ,
  MSD.all.matrix.600.2=foreach(n=1:10,.combine=rbind)%do%{ 
    atom.600.1%>%filter(atom.id==n)%>%select(xu,yu,zu)
  }
  ,
  MSD.all.matrix.600.3=foreach(n=1:10,.combine=rbind)%do%{ 
    atom.600.1.data_table%>%.[atom.id==n]%>%.[,.(xu,yu,zu)]%>%MSD
  }
  ,times = 10
)
timeRecordB(output_message = "benchmark data.table")
## Using data.table gives an 18 times speed-up 


## check if the data.table methods give the same answer as the the dplyr

(  (atom.600.1.data_table%>%.[atom.id==1000]%>%.[,.(xu,yu,zu)]%>%MSD)==(atom.600.1%>%filter(atom.id==1000)%>%select(xu,yu,zu)%>%MSD) )%>%sum()

  



## Next, let us see if the foreach function takes a long time 
microbenchmark(atom.600.1%>%filter(atom.id==100)%>%select(xu,yu,zu)%>%MSD
               ,atom.600.1.data_table%>%.[atom.id==100]%>%.[,.(xu,yu,zu)]%>%MSD
)


## Let us look at the MSD function to see if it takes large portion of time to do. 

