files<-list.files("Documents/PolyembryonyWithYaniv/BurnIn1_0_2_1_0/")

#load("Documents/PolyembryonyWithYaniv/BurnIn-1_-1_0.5_1_0/BurninGenome_w-1_h-1_U0.5_t1_S0_i1")

# summary(z)
# summary(z_burn)
# summary(z2)
# 
# names(z)
# 
# names(z$gen.summary)
# plot(z$gen.summary$gen,z$gen.summary$E_perdiploidgenome)


all.summary <-  do.call(rbind,lapply(files, function(i){
  print(i)
  load(sprintf("Documents/PolyembryonyWithYaniv/BurnIn1_0_2_1_0/%s",i))
  z$gen.summary  <- z$gen.summary %>% mutate(burn = str_remove(i,"Documents/PolyembryonyWithYaniv/BurnIn1_0_2_1_0/BurninGenome_w1_h0_U2_t1_S0_"))
  return(z$gen.summary)
}))

#plot(z$gen.summary$gen,z$gen.summary$E_perdiploidgenome)
ggplot(all.summary, aes(x = gen, y = E_perdiploidgenome, color = burn))+
  geom_line()+
  theme_light()

names(all.summary)

ggplot(all.summary, aes(x = gen, y = mean_w_late_survivors,color = burn))+
  geom_line(alpha = .2)+
  geom_smooth(method = "loess")+
  theme_minimal()+
  ylim(c(.8,.875))

#Burn-in expectation
Neff<-mean(all.summary$two_n[c(1001:2000,3001:4000)])/2


ggplot(all.summary, aes(x = gen, y = E_perdiploidgenome , color = burn))+
  geom_line()+
  theme_minimal() + geom_hline(yintercept=(2/2)*0.5*sqrt(2*pi*500), linetype="dashed", color = "red")

ggplot(all.summary, aes(x = gen, y = two_n , color = burn))+
  geom_line()+
  theme_minimal()


ggplot(all.summary, aes(x = gen, y = cor_w_early_late_all , color = burn))+
  geom_line()+
  theme_minimal()

#Convergence

library("coda")

wpara<-as.mcmc(z$gen.summary[,4:17], start=100, end=2000, thin=50)
g<-geweke.diag(wpara, frac1=0.1, frac2=0.5)
2*pnorm(-abs(g$z))


wpara<-as.mcmc(z$gen.summary[,4:17], start=1000, end=2000)
g<-geweke.diag(wpara, frac1=0.1, frac2=0.5)
2*pnorm(-abs(g$z))


#Tanja's old fashined ways


for (i in files){
  print(i)
}
load("Documents/PolyembryonyWithYaniv/BurnIn-1_-1_0.5_1_0/BurninGenome_w-1_h-1_U0.5_t1_S0_i1")
plot(z$gen.summary$gen,z$gen.summary$E_perdiploidgenome)

plot(z$gen.summary$gen,z$gen.summary$mean_w_late_survivors)
plot(z$gen.summary$gen,z$gen.summary$mean_w_early_survivors)

