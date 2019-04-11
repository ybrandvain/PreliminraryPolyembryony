

### How to handle output
## Summaries
freqP <- function(output){
  # Find the fequency of polyembryony allele over generations
  output$gen.summary                            %>% 
    tidyr::unnest(muts, .drop = TRUE)           %>%
    dplyr::filter(timing == "D")                %>%  
    tidyr::spread(key = s, value = n, fill = 0) %>%
    dplyr::group_by(gen)                        %>%
    dplyr::summarise(p_poly = `11` / (`10` + `11`))
}

mutfreq <- function(output){
  output$gen.summary                            %>% 
  tidyr::unnest(muts, .drop = TRUE)             %>%
  filter(!s %in% c(10,11))
}

numUnique <- function(output){
  # Find the fequency of polyembryony allele over generations
  output$gen.summary                            %>% 
    tidyr::unnest(muts, .drop = TRUE)           %>%
    filter(!s %in% c(10,11))                    %>%
    dplyr::group_by(gen)                        %>%
    dplyr::summarise(n_uqiue_muts = n() )
}



selfingRealized <- function(output){
  primary.self <-  output$gen.summary      %>% 
    tidyr::unnest(w, .drop = TRUE)         %>% 
    dplyr::group_by(gen)                   %>%
    dplyr::summarise(primary_selfing =mean(self))
  realized.self <-  output$gen.summary      %>% 
    tidyr::unnest(w, .drop = TRUE)         %>% 
    dplyr::filter(chosen)                  %>%
    dplyr::group_by(gen)                   %>%
    dplyr::summarise(realized_selfing =mean(self))
  left_join(primary.self, realized.self, by = "gen")
}



fitnessTime <- function(output, focal = "all"){
  # this code is so ugly.. don't look at it
  output.gen.summary <-   output$gen.summary      %>% 
    tidyr::unnest(w, .drop = TRUE)      
  if (focal=="mono") {output.gen.summary <-  output.gen.summary  %>% filter(.,mono == 2) }
  if (focal=="poly") {output.gen.summary <-  output.gen.summary  %>% filter(.,mono == 0) }
  recover()
  # all
  w.summary <- output.gen.summary %>% 
    group_by(gen)    %>%
    summarize(mean.early = mean(w_early), mean.late  = mean(w_late)) %>% left_join(.,
  # chosen
  output.gen.summary %>% 
    filter(chosen)   %>%
    group_by(gen)    %>%
    summarize(mean.early.chosen = mean(w_early), mean.late.chosen  = mean(w_late)), by = "gen" ) %>% left_join(.,
  # not chosen
  output.gen.summary %>% 
    filter(!chosen)   %>%
    group_by(gen)    %>%
    summarize(mean.early.notchosen = mean(w_early), mean.late.notchosen  = mean(w_late)), by = "gen" ) %>%  left_join(.,
  # self
  output.gen.summary %>% 
    filter(self)   %>%
    group_by(gen)    %>%
    summarize(mean.early.self = mean(w_early), mean.late.self  = mean(w_late)) , by = "gen" ) %>%  left_join(.,
  # out
  output.gen.summary %>% 
    filter(!self)   %>%
    group_by(gen)    %>%
    summarize(mean.early.out = mean(w_early), mean.late.out  = mean(w_late)) , by = "gen" ) %>%  left_join(.,
  # self chosen
  output.gen.summary %>% 
    group_by(gen)    %>%
    filter(self & chosen)   %>%
    summarize(mean.early.selfchosen = mean(w_early), mean.late.selfchosen  = mean(w_late)), by = "gen" ) %>%  left_join(.,
  # self notchosen
  output.gen.summary %>% 
    filter(self & !chosen)   %>%
    group_by(gen)    %>%
    summarize(mean.early.selfnotchosen = mean(w_early), mean.late.selfnotchosen  = mean(w_late)) , by = "gen" ) %>%  left_join(.,
  # out chosen
  output.gen.summary %>% 
    filter(!self & chosen)   %>%
    group_by(gen)    %>%
    summarize(mean.early.outchosen = mean(w_early), mean.late.outchosen  = mean(w_late)) , by = "gen" )%>%  left_join(.,
  # out notchosen
  output.gen.summary %>% 
    filter(!self & !chosen)   %>%
    group_by(gen)    %>%
    summarize(mean.early.outnotchosen = mean(w_early), mean.late.outnotchosen  = mean(w_late)), by = "gen" ) 
  w.summary 
}





## Just give me the fucking data
fitnessDataFrame <- function(output){ tidyr::unnest(output$gen.summary, w, .drop = TRUE) %>% data.frame()}
genosDataFrame   <- function(output){ tidyr::unnest(output$gen.summary, muts, .drop = TRUE) %>% data.frame()}
fitnessTime(output = z, focal = "all")
fitnessTime(output = z, focal = "mono")
fitnessTime(output = z, focal = "poly")
freqP(z)
mutfreq(z) 
numUnique(z)
# selfingRealized(z)
# fitnessDataFrame(z) 
# genosDataFrame(z)
