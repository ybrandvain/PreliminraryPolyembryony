library(tidyr)
library("coda")

setwd("Documents/PolyembryonyWithYaniv/BurnIn-1_-1_0.5_1_0")

files<-list.files()
class(files)
geweke_pvalues<-do.call(rbind,lapply(files[1:4], function(i){
  print(i)
  load(i)
  wpara<-as.mcmc(z$gen.summary[,c(4:9,13:14)], start=1000, end=2000, thin=10)
  g<-geweke.diag(wpara, frac1=0.1, frac2=0.5)
  p<-2*pnorm(-abs(g$z))
  pars<-z$params
  return(c(p,as.vector(pars)))
}))

#Visualization

geweke_df<-as.data.frame(geweke_pvalues)
geweke_long <- gather(geweke_df, parameter, p_value, realized_selfing:mean_w_late_all, factor_key = TRUE)
geweke_long = geweke_long %>% unnest()
p<-ggplot(geweke_long, aes(x=parameter, y=p_value)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("GewekePvalues.pdf")

p

dev.off()

write.csv(geweke_long,"GewekePvalues.csv")
