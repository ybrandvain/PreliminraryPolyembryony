# Simulating the evolution of polyembryony
# March 19 2019
# Yaniv Brandvain for collaborative project with Tanja P & Alex Harkness
library(tidyverse)

##########################
# How to live 
##########################
initializeGenomes <- function(n.inds, genomes = NULL){
  if(!is.null(genomes)){return(genomes)}
  tibble(ind    = rep(1:n.inds,2), 
         id     = 10,    # 10 means monoembryony. 11 means poly
         s      =  0,
         h      =  1,    # might change for convinience.. currently meaningless
         timing = "D")
}
introducePoly     <- function(genomes, polyemb.p0){
  genomes %>% 
    dplyr::mutate(id = case_when(id %in% 10:11 ~ 
                                   sample(x = as.double(10:11),
                                          size = n(),
                                          replace = TRUE,
                                          prob = c(1 - polyemb.p0, polyemb.p0)), #10 means monoembryony. 11 means poly
                                 !id %in% 10:11 ~ id))
}
addMutations      <- function(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds){
  getVal <- function(thing, num, this.min = 0, this.max = 1, prelim.vals = NULL){
    if(is.null(prelim.vals)){ prelim.vals <- runif(n = num, min = this.min, this.max) }
    if(thing == "uniform") {   return(prelim.vals)  }
    if(is.numeric(thing))  {   return(rep(thing, num))}
    recover() # we can transform normal to an dist here
  }
  inds   <-   unique(tmp.genomes$ind)
  n.muts <-   rpois(1, U * length(inds) )
  bind_rows(tmp.genomes,                                                     # old genome
            tibble(ind    = sample(inds, size = n.muts, replace = TRUE),     # assigning muts to inds
                   timing = sample(x       = names(dist.timing),             # timing of selection
                                   size    = n.muts,                         # number of muts defined
                                   replace = TRUE,                           # obviously
                                   prob    = dist.timing),                   # right now all muts equi-probable. can change this
                   id = getVal(thing = "uniform", num = n.muts, this.min = 4/n.inds),
                   s  = getVal(thing = fitness.effects, num = n.muts, this.min = 4/n.inds, prelim.vals = id),
                   h  = getVal(thing = dom.effects, num = n.muts))  %>%    # s from uniform as described in ms. can change
              dplyr::mutate(s = ifelse(timing == "B", (1-sqrt(1-s)),s)))   
}
getFitness        <- function(tmp.genomes, dev.to.exclude, adult = TRUE){
  ind.genomes <- tmp.genomes  %>% ungroup()                           %>%
    dplyr::mutate(dup = as.numeric(duplicated(tmp.genomes)))          %>%
    dplyr::filter(!duplicated(tmp.genomes, fromLast = TRUE) & 
                    !timing %in% dev.to.exclude)                               %>%
    dplyr::mutate(w.loc = 1 - (1-dup) * h * s - dup *s)       
  if(adult) {ind.genomes  <- ind.genomes %>% mutate(mono = NA)%>% group_by(ind)}
  if(!adult){ind.genomes  <- ind.genomes %>% group_by(mating, embryo)}
  ind.genomes                                                         %>%   
    dplyr::summarise(w = prod(w.loc), mono = mean(mono))    
}
findMates         <- function(adult.fitness, selfing.rate, n.inds, epsilon = 1e-16){
  outbred.parents <- replicate(3,with(adult.fitness, sample(ind, size = n.inds, replace = TRUE, prob = (w + epsilon)))) 
  self <- data.frame(matrix(rbinom(n = n.inds * 2, size = 1, prob = selfing.rate),ncol = 2))
  colnames(outbred.parents) <- c("mom","dad1","dad2")
  as_tibble(outbred.parents) %>%
    mutate(dad1 = ifelse(self$X1 == 1,mom,dad1), # selfing
           dad2 = ifelse(self$X2 == 1,mom,dad2), # selfing
           mating = 1:n.inds)                    # indexing
}
embryoFitness     <- function(tmp.embryos){
  tmp.embryos                                                 %>% 
    dplyr::group_by(mating)                                   %>%
    dplyr::mutate(mono = sum(parent == "mat" & id == 10)/2)   %>%       
    #dplyr::filter( !(mono == 1 & embryo == "e2") )           %>% 
    dplyr::group_by(embryo, add = TRUE)                       %>% ungroup() %>%
    select(- parent)                                          %>%
    getFitness(dev.to.exclude = "L", adult = FALSE)           %>% ungroup() 
}
favoriteChild     <- function(temp.kidsW, equalizedW = TRUE, compete = TRUE, epsilon = 1e-16){
  temp.kidsW <- temp.kidsW                                    %>%
    dplyr::mutate(alive = rbinom(n = n(),size = 1,prob = w ) )
  if(equalizedW){
    temp.kidsW <- temp.kidsW   %>% 
      group_by(mating) %>%
      mutate(alive = ifelse(mono == 1 , max(alive),  alive))     %>%
      ungroup()
  }
  temp.kidsW <- temp.kidsW                                       %>% 
    dplyr::filter(alive == 1)                                    %>% 
    dplyr::filter( !( mono == 1 & embryo == "e2") )              %>%
    dplyr::group_by(mating) 
  if(!compete){
    temp.kidsW <- temp.kidsW  %>%
      mutate(w = max(w))
  }
  temp.kidsW   %>% 
    sample_n(1,weight = w + epsilon )                         %>% ungroup()  
}
grabInds          <- function(selectedEmbryos, embryos){
  embryoId  <- dplyr::mutate(selectedEmbryos, winners = paste(mating,embryo)) %>% 
    dplyr::select(winners) %>% 
    pull()
  embryos %>% 
    dplyr::mutate(combo = paste(mating, embryo)) %>% 
    dplyr::filter(combo %in% embryoId) %>%
    dplyr::select(ind = mating, id = id, s = s, h = h, timing = timing) 
}
doMeiosis         <- function(tmp.genomes, parents){
  to.meios <- nest(tmp.genomes, -ind)  %>% 
    dplyr::slice(parents)                     %>% 
    dplyr::mutate(mating = 1:length(parents)) %>%
    dplyr::select(mating, data)               %>%
    unnest()           
  poly.allele     <- to.meios$id %in% 10:11 
  to.transmit.hom <- duplicated(to.meios, fromLast = FALSE) & !poly.allele
  pick.one        <- duplicated(to.meios, fromLast = TRUE) | poly.allele  | to.transmit.hom
  # first tranmit het alleels at random
  random.alleles <- to.meios                                    %>% 
    dplyr::filter(!pick.one)                                           %>%
    dplyr::filter(rbinom(n = sum(!pick.one), size = 1, prob = .5) == 1)
  # next be sure to transmit exactly one hom allele
  hom.transmit   <-  to.meios                                   %>% 
    dplyr::filter(to.transmit.hom)
  # finally, pick one of the alelles at the developement locus
  dev.transmit <- to.meios                                      %>% 
    dplyr::filter(poly.allele)                                         %>%
    dplyr::group_by(mating)                                            %>%
    sample_n(size = 1, replace = FALSE)                         %>%
    ungroup()
  # now shove them all together
  bind_rows(random.alleles, hom.transmit, dev.transmit )        %>% 
    nest(-mating)                                               %>% 
    dplyr::arrange(mating) 
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
summarizeGen      <- function(tmp.genomes, mates, embryos, selectedEmbryos){
  selected  <- selectedEmbryos %>% mutate(p = paste(mating, embryo)) %>% select(p)%>% pull 
  muts      <- tmp.genomes %>% mutate(s = ifelse(id %in%  c(10,11), id,s)) %>% group_by(id,s,h,timing) %>% tally() %>% ungroup()
  w.summary <- tibble(
    early_w = getFitness(tmp.genomes,dev.to.exclude = "E", adult = TRUE) %>% select(w) %>%pull(),
    late_w = getFitness(tmp.genomes,dev.to.exclude = "L", adult = TRUE) %>% select(w) %>%pull()) %>%
    summarise(mean_w_late = mean(late_w), mean_w_early = mean(early_w), 
              cor_w_early_late = cor(early_w,late_w))
  list(genome = tmp.genomes %>% mutate, 
       summaries = bind_cols(nest(muts),nest(w.summary)) %>% 
         select(muts = data, w = data1))
}

##########################
##########################
##########################

# running one generation
oneGen <- function(tmp.genomes, n.inds, selfing.rate, U, fitness.effects, dom.effects, dist.timing, equalizedW, compete, just.return.genomes){
  tmp.genomes     <- addMutations(tmp.genomes, U, fitness.effects, dom.effects, dist.timing, n.inds)
  adult.fitness   <- getFitness(tmp.genomes, dev.to.exclude = "E")                       # Adult Fitness
  mates           <- findMates(adult.fitness, selfing.rate, n.inds = n.inds)             # Mating / Selection
  embryos         <- makeBabies(tmp.genomes, mates)                                      # meiosis is in here too
  kidsW           <- embryoFitness(embryos)
  selectedEmbryos <- favoriteChild(kidsW, equalizedW = equalizedW, compete = compete )                                                # pick your child !
  tmp.genomes     <- grabInds(selectedEmbryos = selectedEmbryos, embryos = embryos) %>%  # extract the genomes of selected embryos from our chosen children
    mutate(ind = as.numeric(factor(rank(ind, ties.method = "min"))))                     # ugh.. this last line is kinda gross. but necessary. in means inds are numbered 1:n... this is importnat for  other bits above
  if(just.return.genomes){return(list(genome = tmp.genomes))}
  summarizeGen(tmp.genomes, mates, embryos, selectedEmbryos)
}
# running for a bunch of generations
runSim <- function(n.inds = 1000, selfing.rate = 0, U = .5, fitness.effects  = "uniform", 
                   dom.effects = "uniform", n.gen  = 1000, dist.timing  = c(E = 1/2, B = 0, L = 1/2), 
                   equalizedW = TRUE, compete = TRUE ,
                   introduce.polyem = Inf, polyemb.p0  = .01, genomes = NULL, genome.id = NA,
                   gen.after.loss   = 1,gen.after.fix    = 1, just.return.genomes = FALSE){
  # n.inds           =      1000, 
  # selfing.rate     =         0, # recall selfing = 0 is RANDOM MATING and DOES NOT PRECLUDE SELFING
  # U                =         1, 
  # fitness.effects  = "uniform" or -1, mean uniform, any other number is a fixed value  # need to implement more options. currently takes a fixed val or "uniform"
  # dom.effects      = "uniform",or -1, mean uniform # need to implement more options. currently takes a fixed val or "uniform"
  # n.gen            =      1000, # prob should add an option like "until lost/fixed"
  # dist.timing      = c(E = 1/3, B = 1/3, L = 1/3), add what you like! 
  #                  = numeric options: 1 = c(E = 1/2, B = 0, L = 1/2)
  #                                     2 = c(E = 1, B = 0, L = 0) 
  #                                     3 = c(E = 0, B = 0, L = 1)
  # equalizedW       = TRU. Eshould the expected number of embryos produced by mono and poplyembryonic genos be equivalent? Achieved by group sel at level of mom
  # compete          = compete = TRUE , should embryos be chosen at random or with respect to their fitnesses? 
  # introduce.polyem = Inf       , # gen at which we introduce polyembryony allele
  # polyemb.p0       = .01       , # freq of polyembryony allele once introduced
  # genomes          = NULL        # An option to hand genomes from a previous run
  # genome.id        = NULL 
  # gen.after.loss   = 1
  # gen.after.fix    = 1
  if(fitness.effects == -1){fitness.effects <- "uniform"}
  if(dom.effects == -1){dom.effects <- "uniform"}
  if(length(dist.timing) == 1){
    dist.timing <- list(c(E = 1/2, B = 0, L = 1/2), c(E = 1, B = 0, L = 0), c(E = 0, B = 0, L = 1))[[dist.timing]]
  }
  g             <- 0
  g.after.fix   <- 0 
  g.after.loss  <- 0 
  #g.since.fixed <- 0
  ans           <- list(genome = initializeGenomes(n.inds, genomes)) # Make genomes   # will need to keep track of things... but what?
  keep.going = TRUE
  gen.summary <- list()  
  while(keep.going){  # or stopping rule tbd   # i realize this should be a for loop, but sense that a while loop will give me flexibility for broader stopping rules
    if(g == introduce.polyem){ans$genome <- introducePoly(ans$genome, polyemb.p0)} # introduce polyembryony allele
    g                 <- g + 1
    ans               <- oneGen(ans$genome, n.inds, selfing.rate, U, fitness.effects, 
                                dom.effects, dist.timing, equalizedW = equalizedW, compete = compete, just.return.genomes = just.return.genomes)
    gen.summary[[g]]  <- ans$summaries
    #print(g)
    status <- t(ans$genome  %>% filter(timing == "D")  %>% summarise(loss = sum(id == 11) == 0 , fix = sum(id == 10) == 0) )
    g.after.fix    <- g.after.fix   + as.numeric(status["fix",])
    g.after.loss   <- g.after.loss  + as.numeric(status["loss",])
    if(g >= n.gen   &  (g.after.loss >=  gen.after.loss)   |   (g.after.fix >=  gen.after.fix)){keep.going = FALSE} 
  }
  
  fixed <- ifelse(introduce.polyem == Inf, NA, ifelse(g.after.fix >0, TRUE, FALSE))
  gen.after.fixed.or.lost <- ifelse(introduce.polyem == Inf, NA, ifelse(g.after.fix >0, g.after.fix, g.after.loss ))

  params <- data.frame(n.inds = n.inds, selfing.rate = selfing.rate, U = U, fitness.effects = fitness.effects,
    dom.effects = dom.effects, n.gen = n.gen, g = g, 
    dist.timing = paste(round(dist.timing, digits = 2), collapse = ":"),
    introduce.polyem = introduce.polyem, polyemb.p0  = polyemb.p0 , 
    existing.genome = !is.null(genomes), genom.id = genome.id,  last.gen = g, 
    gen.after.fixed.or.lost  = gen.after.fixed.or.lost, fixed = fixed, equalizedW = equalizedW, compete = compete)
  
  if(just.return.genomes){  return( list( genome = ans$genome, params = params )) }
  gen.summary <- do.call(rbind, gen.summary) %>% mutate(gen = 1:g)
  return(list(genome = ans$genome, gen.summary = gen.summary,params = params))
}
# z <-runSim(n.gen = 1, fitness.effects = 1, dom.effects = -1 ,  gen.after.loss = 15,  gen.after.fix = 15 , polyemb.p0 = 0, introduce.polyem = Inf, just.return.genomes = FALSE)



















olDsummarizeGen      <- function(tmp.genomes, mates, embryos, selectedEmbryos){
  selected  <- selectedEmbryos %>% mutate(p = paste(mating, embryo)) %>% select(p)%>% pull 
  muts      <- tmp.genomes %>% mutate(s = ifelse(id %in%  c(10,11), id,s)) %>% group_by(id,s,h,timing) %>% tally() %>% ungroup()
  w.summary <- left_join(
    embryos                              %>%
      dplyr::group_by(mating)                                 %>%
      dplyr::mutate(mono = sum(parent == "mat" & id == 10)/2)    %>%
      getFitness(dev.to.exclude = "L", adult = FALSE)  %>% ungroup() %>%
      mutate(w_early = w) %>% select(-w)       ,
    embryos                                            %>%
      dplyr::group_by(mating)                                 %>%
      dplyr::mutate(mono = sum(parent == "mat" & id == 10)/2)    %>% ungroup() %>%
      getFitness(dev.to.exclude = "E", adult = FALSE)  %>% ungroup() %>%
      dplyr::mutate(w_late = w) %>% select(-w), 
    by = c("mating", "embryo", "mono"))               %>% 
    left_join(
      mates %>%
        dplyr::mutate(e1 = mom == dad1, e2 = mom == dad2)     %>%
        dplyr::select(- mom , -dad1 , -dad2)                  %>% 
        gather(key = embryo, value = self, -mating),
      by = c("mating", "embryo"))                      %>% 
    left_join(
      embryos                                              %>%
        dplyr::group_by(mating, embryo)                     %>%
        summarise(n_E = sum(timing == "E")  ,
                  n_B = sum(timing == "B") , 
                  n_L = sum(timing == "L") ) ,       
      by = c("mating", "embryo"))         %>%
    dplyr::mutate(z = paste(mating, embryo), 
                  chosen = z %in% selected )                  %>%
    dplyr::select(-z)      
  list(genome = tmp.genomes %>% mutate, 
       summaries = bind_cols(nest(muts),nest(w.summary)) %>% 
         select(muts = data, w = data1))
}