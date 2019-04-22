# run a bunch of invasions (we might want to mess with this more)

simID <-  "BurninGenome_RecessiveLethalEarlyLate1"
load(simID)
#genome <- list(genome = genome, params = z$params) # we wont do this usually, as we will save params  with output


intro.poly.run <- function(simID, n.introductions, n.gen.after.fix, n.times.to.track.fix = 10, n.initial.poly.copies = 1){
  #load(simID)
  poly.models       <- data.frame( equalizedW = c(T,T,F,F),  compete = c(T,F,T,F))
  polyemb.p0        <- as.numeric(n.initial.poly.copies / (2 * genome$params["n.inds"]))
  fitness.effects   <- genome$params["fitness.effects"][1,1]
  dom.effects       <- genome$params["dom.effects"][1,1]
  poly.models.list <- list(equalizedW_compete = 1, equalizedW_nocompete = 2, nonequalizedW_compete = 3,nonequalizedW_nocompete = 4)
  fixed.details <- lapply(poly.models.list,function(X){
    print(X)
    poly.models[X,]
    fixed.outcome <- list()
    times.fixed   <- 0
    sim.summary   <- data.frame(matrix(nrow = n.introductions, ncol = length(genome$params)))
    colnames(sim.summary) <- colnames(genome$params)
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
           just.return.genomes = times.fixed >= n.times.to.track.fix)
    sim.summary[i,] <- this.sim$params
    if(this.sim$params[1,"fixed"]){     
      times.fixed <- times.fixed + 1
      fixed.outcome[[times.fixed]] <- this.sim
      }
    }
    return(list(sim.summary = sim.summary,fixed.outcome = fixed.outcome))
  })
  sim.summary   <- as_tibble(do.call(rbind,lapply(intro.results, function(X){X$sim.summary})))
  fixed.details <- bind_rows(lapply(fixed.details, function(TMp){
    bind_rows(lapply(TMp$fixed.outcome, function(TMP){
      bind_cols(as_tibble(nest(TMP$gen.summary)),
                as_tibble(nest(TMP$genome)),
                as_tibble(TMP$params))%>% 
        rename( gen.summary = data , genome = data1)
    }))
  }))
  return(list(sim.summary = sim.summary, fixed.details = fixed.details))
}


intro.poly.run(simID = "BurninGenome_RecessiveLethalEarlyLate1", n.introductions = 10, n.gen.after.fix = 500)
