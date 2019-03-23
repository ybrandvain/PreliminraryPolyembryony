

### How to handle output
## Summaries
freqP <- function(output){
  # Find the fequency of polyembryony allele over generations
  recover()
  output$gen.summary                            %>% 
    tidyr::unnest(muts, .drop = TRUE)           %>%
    dplyr::filter(timing == "D")                %>%  
    tidyr::spread(key = s, value = n, fill = 0) %>%
    dplyr::group_by(gen)                        %>%
    dplyr::summarise(p_poly = `11` / (`10` + `11`))
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






## Just give me the fucking data
fitnessDataFrame <- function(output){ tidyr::unnest(output$gen.summary, w, .drop = TRUE) %>% data.frame()}
genosDataFrame   <- function(output){ tidyr::unnest(output$gen.summary, muts, .drop = TRUE) %>% data.frame()}

fitnessTime(z)
# freqP(z)
# selfingRealized(z)
# fitnessDataFrame(z) 
# genosDataFrame(z)
