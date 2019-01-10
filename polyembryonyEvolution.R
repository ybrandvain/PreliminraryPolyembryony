# Yaniv Brandvain & Tanja Pyhajarvi
# The evolution of polyembyony 
# Last updated 1/19


#########################################################################################
# Required packages                                                                     #
library(ggplot2); library(forcats); library(dplyr)                                      #
# Helpful for making figures                                                            #
gg_color_hue <- function(n) {                                                           #
  hues = seq(15, 375, length = n + 1)                                                   #
  hcl(h = hues, l = 65, c = 100)[1:n]                                                   #
}                                                                                       #
#########################################################################################

# Background: Some plants [including pines] have > 1 embryo / seed, but only one embryo survives
#             Here, we explore the evolution of "Polyembryony"
#             First, we ask if polyembryony can evolve as a mechanism to prevent selfing.
#               We assume a finit population size and an infinite sites model for deleterious recessives
#               Fitness is multiplicative. The input DFE can be customized... but is currently either a fixed s or a uniform
#               We assume that selection acts equivavlently on loci during embryonic and early zygotic selection


#######################################################
#Primary Functions                                    #



# Calcualting individual fitness 
indInfo <- function(my.inds, s, time.sel = "both"){
  # find fitness etc of each ind
  # currently assumes multiplicative fitness with a fixed selective coefficient, and full recessivity.
  indFitness <- function(IND, s){
    if(length(IND) == 0){return(c(w = 1, n.mut = 0))}
    homs    <- duplicated(IND)
    n.hom   <- sum(homs)
    n.mut   <- sum(!homs)
    if(s == "uniform"){
      # Here i try to select by time of actions..  THIS NEEDS WORK / THOUGHT
      early   <- IND < 0  
      if(time.sel == "both") {sel.efx <- abs(IND[homs])}           # Assumes all loci selected between embryos and between seed(lings)
      if(time.sel == "early"){sel.efx <- abs(IND[homs & early]) }  # Looks at loci only selected between embryos
      if(time.sel == "late") {sel.efx <- abs(IND[homs & !early]) } # Looks at loci only selected between seed(lings)
      fitness <- prod( 1 - sel.efx/1000 ) 
      if(is.na(fitness)){fitness <- 1}
    } # here fitness comes from the random uniform [given arbitrary id]
    if(is.numeric(s)){ fitness <- (1 - s)^n.hom } # here we have a fixed s
    return(c(w = fitness, n.mut = n.mut))
  }
  all.fitness <- sapply(my.inds, indFitness, s = s)
  return(all.fitness)
}




# Mating and selection 
selectANDMate <- function(tmp.inds, w, n.inds, selfing.rate = "random"){
  # Adds selfing rate on otp of incidental selfing... so no explicit avoidance of self matings
  # Makes two embryos / seed. "Polyembryony"  mean we select between these embryos, rather than picking one at random
  ind.id <- seq_along(tmp.inds)
  tmp.chroms <- data.frame(tmp.mat = sample(ind.id, replace = TRUE, prob = (w + 1e-10)),  # so we're not all dead
                           tmp.pat1 = sample(ind.id, replace = TRUE, prob = (w + 1e-10)),  #  so we're not all dead
                           tmp.pat2 = sample(ind.id, replace = TRUE, prob = (w + 1e-10)))  #  so we're not all dead
  if(is.numeric(selfing.rate)){
    # currently ignores the nonindependence of dads.. can induce an autocorrelation if we wish..
    # also currently models polyembryony from same haploid mat... can change / compare [with some effort]
    selfed1 <- sample(c(TRUE,FALSE), size = n.inds, replace = TRUE, prob = c(selfing.rate, 1-selfing.rate))
    tmp.chroms$tmp.pat1[selfed1] <- tmp.chroms$tmp.mat[selfed1]
    selfed2 <- sample(c(TRUE,FALSE), size = n.inds, replace = TRUE, prob = c(selfing.rate, 1-selfing.rate))
    tmp.chroms$tmp.pat2[selfed2] <- tmp.chroms$tmp.mat[selfed2]
  }
  return(tmp.chroms)
}


# Synagamy and meiosis (for both embryos)
syngamyAndMeiosis <- function(mat, pat1, pat2, n.muts){
  doMeiosis <- function(chrom){
    homs <- duplicated(chrom)
    hom.alleles <- chrom[homs]
    het.alleles <- setdiff(chrom, hom.alleles)
    to.keep     <- sample(c(TRUE, FALSE), length(het.alleles), replace = TRUE)
    transmitted <- c(het.alleles[to.keep], hom.alleles)
    return(transmitted)
  }
  mat.gamete <- doMeiosis(mat)
  pat1.gamete <- doMeiosis(pat1)
  pat2.gamete <- doMeiosis(pat2)
  new.muts   <- round(runif(n.muts, min = -1000, max = 1000), digits = 10)
  # because only one kid will survive and bc all mutations are unique im just develeoping one set of new mutations
  zygote1     <- unlist(c(mat.gamete, pat1.gamete, new.muts))
  zygote2     <- unlist(c(mat.gamete, pat2.gamete, new.muts))
  return(list(zygote1, zygote2))
}




#one gen
oneGen <- function(n.inds, U, selfing.rate = 0, s, tmp.inds, w.vector, poly.geno = NULL, poly.sel.all = FALSE){
  # Here we make a generations. 
  # Population is of size n.inds. 
  # Selfing rate is the probability of selfing above and beyound random mating 
  # s is the recessive fitness cost of new mutations. can be a number or come from a distribution [currently uniform only]
  # Assumes we already has already been assigned fitness as a seed(ling)
  # We start by selecting individuals and finding mates
  mates         <-  selectANDMate(tmp.inds     = tmp.inds,
                                  w            =  w.vector,
                                  n.inds       = n.inds,
                                  selfing.rate = selfing.rate)
  # We then layer on mutations
  muts          <- rpois(n = n.inds, lambda = U)
  # Here everyone is making two embryos.. I decide they are selectively or randomly transmitted below. This is slow and wasteful
  embryos <- mapply(syngamyAndMeiosis, 
                    mat  = tmp.inds[mates$tmp.mat] , 
                    pat1 = tmp.inds[mates$tmp.pat1] , 
                    pat2 = tmp.inds[mates$tmp.pat2] , 
                    n.muts = muts)
  # these ifs are all gross and prolly better treated as function calls ... will deal with later
  if(is.null(poly.geno) & !poly.sel.all){chosen.embryo <- 0 } # functionally, no polyembryony
  if(!is.null(poly.geno) | poly.sel.all){
    # if either all are polyembryonic or polyembryony is evolving
    embryo.info      <- indInfo(my.inds = embryos, s = 1)
    embryo.fitnesses <- matrix(embryo.info["w",],ncol = 2, byrow = TRUE) + 1e-10
    embryo2.fitness  <- embryo.fitnesses[,2] / rowSums(embryo.fitnesses) # here we can allow polyembryony to vary by individual
    if(!is.null(poly.geno)){ 
      # if polyembryony is evolving
      embryo2.fitness <- ifelse(poly.geno[mates$tmp.mat] == 0, .5, embryo2.fitness) 
      # right now assumes dominant effect of maternal diploid . can change
    }
    chosen.embryo    <- rbinom(n = n.inds, size = 1, prob = embryo2.fitness) 
  }
  new.zygotes      <- seq(1,2*n.inds-1,2) + chosen.embryo
  new.inds <- embryos[new.zygotes]
  if(is.null(poly.geno)){return(new.inds)} # return new genos at selected loci only if polyembryony is not evolving
  
  # If polyembryony is evolving we need to find genotypes at polyembryony locus and return these too
  dad <- ifelse(chosen.embryo == 0, mates$tmp.pat1,mates$tmp.pat2) # find the dad
  dad.poly.geno <- poly.geno[dad]
  mom.poly.geno <- poly.geno[mates$tmp.mat]
  new.poly.geno <- (rbinom(n.inds, 1, mom.poly.geno ) + rbinom(n.inds, 1, dad.poly.geno ) ) / 2 # synagmy and meiosis. 
  return(list(new.inds = new.inds, new.poly.geno = new.poly.geno))
}


### Looping across generations 
polyembryonyEvoluationSim <- function(n.inds, gen.intoroduced, U, selfing.rate = 0, s = 1, return.genos = FALSE, 
                                      pe.intro.freq = (1 / (2*n.inds)), pe.intro.tries = 0){
  n.chroms  <- 2 * n.inds
  gen.sum	<- data.frame(mean.fitness = NA, mean.muts = NA,n.unique.muts = NA, 
                        pe.allele.freq = NA, cov.pe.w = NA, cov.pe.nmuts = NA)
  inds      <- replicate(n.inds,list())  
  for (g in 1:gen.intoroduced){
    ind.info      <- indInfo(my.inds = inds, s = s) 
    # There is some redundancy in that im essentially doing this twice.. obvious place for future speed up
    gen.sum[g,1:3]   <- c(rowMeans(ind.info), length(unique(do.call(c,inds))))
    inds          <- oneGen(n.inds, U, selfing.rate, s, tmp.inds = inds,  w.vector = ind.info["w",],poly.geno = NULL)
  }
  inds.before.pe.mut <- inds
  intro.attempt <- 0    # I'm working on allowing more than one introduction per simulation.. not implemented at the moment
  while(intro.attempt < pe.intro.tries){
    poly.geno <- rbinom(n.inds,2, pe.intro.freq)/2
    while(  (mean(poly.geno) * (1-mean(poly.geno))) !=0 ){
      # need to fix generation counting for more than one trial
      g <- g + 1
      ind.info  <- indInfo(my.inds = inds, s = s) 
      # There is some redundancy in that im essentially doing this twice.. obvious place for future speed up
      gen.sum[g,1:3]   <- c(rowMeans(ind.info), length(unique(do.call(c,inds))))
      gen.sum[g,"pe.allele.freq"] <- mean(poly.geno) 
      gen.sum[g,"cov.pe.w"] <- cov(poly.geno, ind.info["w",]) 
      gen.sum[g,"cov.pe.nmuts"] <- cov(poly.geno, ind.info["n.mut",]) 
      inds          <- oneGen(n.inds, U, selfing.rate, s, tmp.inds = inds,  w.vector = ind.info["w",], poly.geno = poly.geno)
      poly.geno     <- inds[["new.poly.geno"]]
      inds          <- inds[["new.inds"]]
    }
    intro.attempt <- intro.attempt + 1
  }
  if(!return.genos){return(data.frame(n.inds = n.inds, selfing.rate = selfing.rate, U = U, gen = 1:g, gen.sum))}
  if(return.genos){return(list(inds = inds, out = data.frame(n.inds = n.inds, selfing.rate = selfing.rate, U = U, gen = 1:g, gen.sum)))}
}
#######################################################

# Working example
this.sim <- polyembryonyEvoluationSim(n.inds = 1000, gen.intoroduced = 250, U = 1, selfing.rate = 0, s = "uniform", pe.intro.tries = 1, pe.intro.freq = .01)
