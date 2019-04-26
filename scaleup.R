# run a bunch of invasions (we might want to mess with this more)
setwd("../../../../../../projects/polyembryony/")
#runs the scripts
source("polyembryonyEvolution.R")
source("processPoly.R")


intro.poly.run <- function(simID, n.introductions, n.gen.after.fix, n.times.to.track.fix = 10, n.initial.poly.copies = 1,which.model = 1:4){
  load(simID)
  poly.models       <- data.frame( equalizedW = c(T,T,F,F),  compete = c(T,F,T,F))
  polyemb.p0        <- as.numeric(n.initial.poly.copies / (2 * z$params["n.inds"]))
  fitness.effects   <- z$params["fitness.effects"][1,1]
  dom.effects       <- z$params["dom.effects"][1,1]
  poly.models.list  <- list(equalizedW_compete = 1, equalizedW_nocompete = 2, nonequalizedW_compete = 3,nonequalizedW_nocompete = 4)
  poly.models.list  <- poly.models.list[which.model] 
  dist.timing       <- as.numeric(strsplit(as.character(z$params["dist.timing"][1,1]),":")[[1]])
  names(dist.timing)<- c("E","B","L")
  #mygenome          <- z$genome
  intro.results <- lapply(poly.models.list,function(X){
    print(X)
    poly.models[X,]
    fixed.outcome <- list()
    times.fixed   <- 0
    sim.summary   <- data.frame(matrix(nrow = n.introductions, ncol = length(z$params)))
    colnames(sim.summary) <- colnames(z$params)
    for(i in 1:n.introductions) { # should be n.introductions
      print(i)
      this.sim <- runSim(n.gen   = 0, 
                         fitness.effects     = fitness.effects, 
                         dom.effects         = dom.effects ,  
                         gen.after.loss      = 1,  
                         gen.after.fix       = ifelse(times.fixed < n.times.to.track.fix, n.gen.after.fix,1), 
                         polyemb.p0          = polyemb.p0, 
                         introduce.polyem    = 0,
                         equalizedW          = poly.models[X,"equalizedW"] ,
                         compete             = poly.models[X,"compete"],
                         genomes             = z$genome,
                         genome.id           = simID,
                         just.return.genomes = times.fixed >= n.times.to.track.fix,
                         dist.timing         = dist.timing)
      sim.summary[i,] <- this.sim$params
      if(this.sim$params[1,"fixed"]){     
        times.fixed <- times.fixed + 1
        fixed.outcome[[times.fixed]] <- this.sim
      }
    }
    return(list(sim.summary = sim.summary,fixed.outcome = fixed.outcome))
  })
  sim.summary   <- as_tibble(do.call(rbind,lapply(intro.results, function(X){X$sim.summary})))
  fixed.details <- bind_rows(lapply(intro.results$fixed.outcome , function(TMp){
    bind_rows(lapply(TMp$fixed.outcome, function(TMP){
      bind_cols(as_tibble(nest(TMP$gen.summary)),
                as_tibble(nest(TMP$genome)),
                as_tibble(TMP$params))%>% 
        rename( gen.summary = data , genome = data1)
    }))
  }))
  return(list(sim.summary = sim.summary, fixed.details = fixed.details))
}

a <-intro.poly.run(simID = "BurnInGenomes/BurninGenome_RecessiveLethalEarlyLateSelfing_0.9_10", n.introductions = 2, n.gen.after.fix = 1)
#intro.poly.run(simID = "BurninGenome_RecessiveLethalEarlyLate1", n.introductions = 100, n.gen.after.fix = 500)
