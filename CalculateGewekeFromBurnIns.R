library(tidyr)
library(coda)
library(ggplot2)
library(dplyr)

#set the Geweke working directory
setwd("/homeappl/home/pyhajar/akaprojlink/Polyembryony/Geweke")

#Read filenames for each Burn In genome
files<-read.table("GenomePaths3000.txt", as.is=T)
files<-files[,1]


#read genome summary values from genome objects (produced by polyembryonyEvolution.R)
geweke_pvalues<-do.call(rbind,lapply(files, function(i){
  print(i)
  load(i)
  # read as an mcmc object
  wpara<-as.mcmc(z$gen.summary[,c("realized_selfing", "primary_selfing","E_perhaploidgenome", "L_perhaploidgenome", "mean_w_late_survivors", "mean_w_early_survivors",
                                  "mean_w_early_all","mean_w_late_all")], start=length(z$gen.summary$gen)/2, end=length(z$gen.summary$gen), thin=10)
  # calculate Geweke diagnostic statistics
  g<-geweke.diag(wpara, frac1=0.5,frac2=0.5)
  
  #returned z.score converted to p-value
  p<-2*pnorm(-abs(g$z))
  
  pars<-z$params
  return(c(p,as.vector(pars)))
}))


save(file="GewekePvalues3000New.Robj", geweke_pvalues)
