# Simulating the evolution of polyembryony
# March 19 2019
# Yaniv Brandvain for collaborative project with Tanaj & Alex Harkness
library(tidyverse)

##########################
# How to live 
##########################
initializeGenomes <- function(n.inds, genomes){
  if(!is.null(genomes)){return(genomes)}
  tibble(ind    = rep(1:n.inds,2), 
         s      = 10,    #10 means monoembryony. 11 means poly
         h      = 1,    # might change for convicence.. currently meaningless
         timing = "D")
}
introducePoly     <- function(genomes, polyemb.p0){
  genomes %>% 
    mutate(s = case_when(s %in% 10:11 ~ 
                           sample(x = 10:11,
                                  size = n(),
                                  replace = TRUE,
                                  prob = c(1 - polyemb.p0, polyemb.p0)), #10 means monoembryony. 11 means poly
                         !s %in% 10:11 ~ s))
}
addMutations      <- function(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds){
  getVal <- function(thing, num, this.min = 0, this.max = 1){
    if(is.numeric(thing)) {return(rep(thing, num))}
    if(fitness.effects == "uniform"){return(runif(n = n.muts, min = this.min, this.max))}
    recover()
  }
  inds   <-   unique(tmp.genomes$ind)
  n.muts <-   rpois(1, U * length(inds) )
  bind_rows(tmp.genomes,                                                     # old genome
            tibble(ind    = sample(inds, size = n.muts, replace = TRUE),     # assigning muts to inds
                   timing = sample(x       = names(dist.timing),             # timing of selection
                                   size    = n.muts,                         # number of muts defined
                                   replace = TRUE,                           # obviously
                                   prob    = dist.timing),                   # right now all muts equi-probable. can change this
                   s = getVal(thing = fitness.effects, num = n.muts, this.min = 4/n.inds),
                   h = getVal(thing = dom.effects, num = n.muts))  %>%    # s from uniform as described in ms. can change
                   mutate(s = ifelse(timing == "B", (1-sqrt(1-s)),s)))   
}
getFitness        <- function(tmp.genomes, dev.to.exclude, adult = TRUE){
  ind.genomes <- tmp.genomes  %>% ungroup()                           %>%
    mutate(dup = as.numeric(duplicated(tmp.genomes)))                 %>%
    filter(!duplicated(tmp.genomes, fromLast = TRUE) & 
           !timing %in% dev.to.exclude)                               %>%
    mutate(s = ifelse(timing == "D", 0, s),
           w.loc = 1 - (1-dup) * h * s - dup *s,
           mono = ifelse(!adult, mean(mono),NA ))        
  if(adult) {ind.genomes  <- ind.genomes %>% group_by(ind)}
  if(!adult){ind.genomes <- ind.genomes %>% group_by(mating, embryo)}
  ind.genomes                                                         %>%   
    summarise(w = prod(w.loc), mono = mean(mono))  
}
findMates         <- function(adult.fitness, selfing.rate, n.inds){
  outbred.parents <- replicate(3,with(adult.fitness, sample(ind, size = n.inds, replace = TRUE, prob = w ))) 
  self <- data.frame(matrix(rbinom(n = n.inds * 2, size = 1, prob = selfing.rate),ncol = 2))
  colnames(outbred.parents) <- c("mom","dad1","dad2")
  as_tibble(outbred.parents) %>%
    mutate(dad1 = ifelse(self$X1 == 1,mom,dad1), # selfing
           dad2 = ifelse(self$X2 == 1,mom,dad2), # selfing
           mating = 1:n.inds)                    # indexing
}
embryoFitness     <- function(tmp.embryos){
  tmp.embryos                                                 %>% 
    group_by(mating)                                          %>%
    mutate(mono = sum(parent == "mat" & s == 10))             %>%       
    filter( !(mono == 2 & embryo == "e2") )                   %>% 
    group_by(embryo, add = TRUE)                              %>% ungroup() %>%
    getFitness(dev.to.exclude = "L", adult = FALSE)           %>% ungroup() 
}
favoriteChild     <- function(temp.kidsW){
  temp.kidsW                                                  %>%
  mutate(alive = rbinom(n = n(),size = 1,prob = w ) )         %>% 
    filter(alive == 1)                                        %>% # here is a place i could add a cost to mon or ploy by fucking with fitness. for ecample i could ensure that the prob survival is the same !!!
    group_by(mating)                                          %>%
    sample_n(1,weight = w)                                    %>% ungroup()
}
grabInds          <- function(selectedEmbryos, embryos){
  embryoId  <- mutate(selectedEmbryos, winners = paste(mating,embryo)) %>% 
    select(winners) %>% 
    pull()
  embryos %>% 
    mutate(combo = paste(mating, embryo)) %>% 
    filter(combo %in% embryoId) %>%
    select(ind = mating, s = s, h = h, timing = timing) 
}
doMeiosis         <- function(tmp.genomes, parents){
  to.meios <- nest(tmp.genomes, -ind)  %>% 
    slice(parents)                     %>% 
    mutate(mating = 1:length(parents)) %>%
    select(mating, data)               %>%
    unnest()           
  poly.allele     <- to.meios$s %in% 10:11 
  to.transmit.hom <- duplicated(to.meios, fromLast = FALSE) & !poly.allele
  pick.one        <- duplicated(to.meios, fromLast = TRUE) | poly.allele  | to.transmit.hom
  # first tranmit het alleels at random
  random.alleles <- to.meios                                    %>% 
    filter(!pick.one)                                           %>%
    filter(rbinom(n = sum(!pick.one), size = 1, prob = .5) == 1)
  # next be sure to transmit exactly one hom allele
  hom.transmit   <-  to.meios                                   %>% 
    filter(to.transmit.hom)
  # finally, pick one of the alelles at the developement locus
  dev.transmit <- to.meios                                      %>% 
    filter(poly.allele)                                         %>%
    group_by(mating)                                            %>%
    sample_n(size = 1)                                          %>%
    ungroup()
  # now shove them all together
  bind_rows(random.alleles, hom.transmit, dev.transmit )        %>% 
    nest(-mating)                                               %>% 
    arrange(mating)
}
makeBabies        <- function(tmp.genomes, mates){   
  bind_cols(
    doMeiosis(tmp.genomes, mates$mom)  %>% select(mating = mating, mat = data),
    doMeiosis(tmp.genomes, mates$dad1) %>% select(e1 = data),
    doMeiosis(tmp.genomes, mates$dad2) %>% select(e2 = data))                       %>%
    gather(key = embryo, value = pat, -mating, - mat)                               %>%
    gather(key = parent, value = haplo, -mating, - embryo)                          %>% # syngamy
    unnest(haplo)                                                              # not sure about unnesting here.... 
}

##########################
##########################
##########################

# running one generation
oneGen <- function(tmp.genomes, n.inds, selfing.rate, U, fitness.effects, dom.effects, dist.timing){
  tmp.genomes     <- addMutations(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds)
  adult.fitness   <- getFitness(tmp.genomes, dev.to.exclude = "E")                       # Adult Fitness
  mates           <- findMates(adult.fitness, selfing.rate, n.inds = n.inds)             # Mating / Selection
  embryos         <- makeBabies(tmp.genomes, mates)                                      # meiosis is in here too
  kidsW           <- embryoFitness(embryos)
  selectedEmbryos <- favoriteChild(kidsW)                                                # pick your child !
  tmp.genomes     <- grabInds(selectedEmbryos = selectedEmbryos, embryos = embryos) %>%  # extract the genomes of selected embryos from our chosen children
    mutate(ind = as.numeric(factor(sort(rank(ind, ties.method = "min")))))               # ugh.. this last line is kinda gross. but necessary. in means inds are numbered 1:n... this is importnat for  other bits above
  return(tmp.genomes)
}
# running for a bunhc of generations
runSim <- function(n.inds = 1000, selfing.rate = 0, U = 1, fitness.effects  = "uniform", 
                   dom.effects = "uniform", n.gen  = 1000, dist.timing  = c(E = 1/3, B = 1/3, L = 1/3), 
                   introduce.polyem = Inf, polyemb.p0  = .01, genomes = NULL){
  # n.inds           =      1000, 
  # selfing.rate     =         0, # recall selfing =0 is RANDOM MATING and DOES NOT PRECLUDE SELFING
  # U                =         1, 
  # fitness.effects  = "uniform", # need to implement more options. currently takes a fixed val or "uniform"
  # dom.effects      = "uniform", # need to implement more options. currently takes a fixed val or "uniform"
  # n.gen            =      1000, # prob should add an option like "until lost/fixed"
  # dist.timing      = c(E = 1/3, B = 1/3, L = 1/3),
  # introduce.polyem = Inf       , # gen at which we introduce polyembryony allele
  # polyemb.p0       = .01       , # freq of polyembryony allele once introduced
  # genomes          = NULL        # An option to hand genomes from a previous run
  genomes      <- initializeGenomes(n.inds, genomes) # Make genomes   # will need to keep track of things... but what?
  g <- 0
  while(g < n.gen){  # or stopping rule tbd   # i realize this should be a for loop, but sense that a while loop will give me flexibility for broader stopping rules
    if(g == introduce.polyem){genomes <- introducePoly(genomes, polyemb.p0)} # introduce polyembryony allele
    g         <- g + 1
    genomes  <- oneGen(genomes, n.inds, selfing.rate, U, fitness.effects, dom.effects, dist.timing)
    # NOTE: I prolly want to return a list with more stuff in future in addition to genomes
  }
}
runSim()