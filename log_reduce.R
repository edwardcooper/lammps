

source("https://raw.githubusercontent.com/edwardcooper/lammps/master/long_to_wide.R")

## Reduce the size of the scraped log.lammps file. 

log_reduce=function(infile,interval,outfile){
  log=read.table(file=infile)
  log=log%>%long_to_wide()
  reduce_index=seq(from=1,to=nrow(log),by=interval)
  log=log[reduce_index,]
  write.table(x=log,file=outfile,sep=",")
}

# log_reduce(infile ="log320_solomon.csv"
#            ,interval=10
#            ,outfile="test.csv" )
# 
# test_csv=read.table("test.csv",sep=",")
# colnames(test_csv)
