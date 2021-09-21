setwd("Q:/Technical/Python")
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggallin)
library(viridis)
library(meta)

### local invasiveness df
fra_inv_file <- "Q:/Technical/Python/popsims/fra_inv"
fra_dis_file <- "Q:/Technical/Python/popsims/fra_dis"
fra_glob_file <- "Q:/Technical/Python/popsims/fra_globinv"
usa_inv_file <- "Q:/Technical/Python/popsims/usa_inv"
usa_dis_file <- "Q:/Technical/Python/popsims/usa_dis"
usa_glob_file <- "Q:/Technical/Python/popsims/usa_globinv"

# function that plots vaccine vs pooled IRR; shape = strategy, colour = age grp
pl <- function(des,metadataframe) {
  ggplot(metadataframe %>% filter(description == des), 
         group = agegrp, 
         aes(color = agegrp, shape = vaxstrategy)) + 
    geom_hline(aes(yintercept = exp(1)), linetype = 'dashed') + 
    geom_errorbar(aes(x = vaxperiod, ymin = fixed.lo, ymax = fixed.hi), width = 0.01,
                  position = position_dodge(width = 0.5)) +
    geom_point(aes(x = vaxperiod, y = fixed), 
               position = position_dodge(width = 0.5)) + ggtitle(des)+
    labs(y = 'log(Pooled IRR)', x = "", color = "", shape = "") +
    scale_y_continuous(trans = pseudolog10_trans) + theme_bw() + theme_minimal() +
    scale_color_viridis(discrete = T)}

pl2 <- function(agegroup,metadataframe,vaxstrat) {
  ggplot(metadataframe %>% filter(agegrp == agegroup, vaxstrategy == vaxstrat), 
         group = description, 
         aes(color = description)) + 
    geom_hline(aes(yintercept = exp(1)), linetype = 'dashed') + 
    geom_errorbar(aes(x = vaxperiod, ymin = fixed.lo, ymax = fixed.hi), width = 0.01,
                  position = position_dodge(width = 0.5)) +
    geom_point(aes(x = vaxperiod, y = fixed), 
               position = position_dodge(width = 0.5)) + ggtitle(agegroup)+
    labs(y = 'log(Pooled IRR)', x = "", color = "", shape = "") +
    scale_y_continuous(trans = pseudolog10_trans) + theme_bw() + theme_minimal() +
    scale_color_viridis(discrete = T)}

pooled_est <- function(df, des, vaxperiod2, agegroup) {
    #data (read, tidy, and structure as needed)
    data_temp <- df %>% filter(description == des) %>% filter(X == vaxperiod2) %>% filter(agegrp == agegroup)
    data_temp <- data_temp %>% filter(low > 0 & high > 0 & IRR > 0) %>% filter(!is.infinite(IRR))
    #data_temp <- data_temp %>% mutate(error = (log(low) - log(IRR))/(-1.96)) # check if error returns same results
  
    #assign model (typically a nice meta. function from one of several packages such as meta, metafor, or netmeta)
    meta <- metagen(log(IRR), lower=log(low), upper = log(high), studlab = simnum, data = data_temp) #fit generic meta-analysis to an object
    #k <- metagen(log(IRR), error, studlab = simnum, data = data_temp)
    
    #viz (draw a standard forest plot or metaregression plot)
    #pdf(paste(des, "_", vaxperiod2, "_", agegroup, ".pdf", sep = ""), width = 15, height = 10) 
    #forest(meta)
    #dev.off()
    
    return(meta) #grid-based graphics so a bit of work to resize
}

IRRfun <- function(foldernam) {
  setwd(foldernam)
  
  countrynam <- gsub(".*popsims/(.+)_.*", "\\1", foldernam)
  vaxstrategy <- gsub(paste(".*", countrynam,"_(.+).*", sep = ""), "\\1", foldernam)
  
  if (vaxstrategy == 'dis') {
    vaxstrat <- 'dis_inc'
  }
  if (vaxstrategy == 'globinv') {
    vaxstrat <- 'glob_inv'
  }
  if (vaxstrategy == 'inv') {
    vaxstrat <- 'local_inv'
  }
  
  temp <- list.files(pattern="*.csv")
  myfiles2 <- lapply(temp, read.csv)
  myfiles <- myfiles2
  
  simnum <- sub("_.*", "", temp)
  simnum <- as.integer(sub("simulation", "", simnum))
  
  vaxperiod <- gsub("_.*", "" ,gsub(".*_p(.+)_.*", "p\\1", temp))
  
  simparams <- read.csv("Q:/Technical/Python/simparams.csv")
  
  for (i in 1:length(myfiles)) {
    myfiles[[i]]$simnum <- simnum[i]
    myfiles[[i]]$vaxperiod <- vaxperiod[i]
    myfiles[[i]]$country <- countrynam
    myfiles[[i]]$description <- simparams[which(simparams$simnum == simnum[i]),]$description
    myfiles[[i]]$vaxstrat <- vaxstrat
  }

  agegrpstatsfiles <- grepl('agegrpstats', temp, fixed = TRUE)
  agegrpstats_temp <- myfiles[agegrpstatsfiles]
  agegrpstats <<- bind_rows(agegrpstats_temp)
  if (vaxstrat != 'local_inv') {
    agegrpstats_prePCV <- agegrpstats_locinv %>% filter(vaxperiod == 'prevacc')
    agegrpstats <<- bind_rows(agegrpstats_prePCV, agegrpstats)
  }
  
  disincfiles <- grepl('disinc.csv', temp, fixed = TRUE)
  disinc_temp <- myfiles[disincfiles]
  disinc <- bind_rows(disinc_temp)
  if (vaxstrat != 'local_inv') {
    disinc_prePCV <- disinc_locinv %>% filter(vaxperiod == 'prevacc')
    disinc <<- bind_rows(disinc_prePCV, disinc)
  }
  
  outputfiles <- grepl('output', temp, fixed = TRUE)
  output_temp <- myfiles[outputfiles]
  output <- bind_rows(output_temp)
  
  VTserofiles <- grepl('VTsero', temp, fixed = TRUE)
  VTsero_temp <- myfiles[VTserofiles]
  VTsero <<- bind_rows(VTsero_temp)
  
  IRRfiles <- grepl('IRR_allperiods', temp, fixed = TRUE)
  IRR_temp <- myfiles[IRRfiles]
  IRR <- bind_rows(IRR_temp)
  
  carrprevfiles <- grepl('carrprev_eachperiod', temp, fixed = TRUE)
  carrprev_temp <- myfiles[carrprevfiles]
  carrprev <<- bind_rows(carrprev_temp)
  
  disinceachfiles <- grepl('disinc_eachperiod', temp, fixed = TRUE)
  disinceach_temp <- myfiles[disinceachfiles]
  disinceach <- bind_rows(disinceach_temp)
  
  carrierstmaxfiles <- grepl('carrierstmax', temp, fixed = TRUE)
  carrierstmax_temp <- myfiles[carrierstmaxfiles]
  carrierstmax <- bind_rows(carrierstmax_temp)
  if (vaxstrat != 'local_inv') {
    carrierstmax_prePCV <- carrierstmax_locinv %>% filter(vaxperiod == 'prevacc')
    carrierstmax <- bind_rows(carrierstmax_prePCV, carrierstmax)
  }
  
  # simulations missing postPCV30 files (agegrpstats and disinc)
  post30_disinc_fin <- unique((disinc %>% filter(vaxperiod == 'postPCV30'))$simnum)
  post30_output_fin <- unique((output %>% filter(vaxperiod == 'postPCV30'))$simnum)
  post30missing <- setdiff(post30_output_fin, post30_disinc_fin)
  
  for (i in post30missing) {
    output_t <- output %>% filter(simnum == i, vaxperiod == 'postPCV30') %>%
      filter(X != 'time') %>% filter(X != 'sim_time')
    colnames(output_t)[1:2] <- c('Serotype', 'carrprev')
    indivs <- 5000
    output_t$carrprev <- output_t$carrprev/indivs
    carriers_t <- carrierstmax %>% filter(simnum == i, vaxperiod == 'postPCV30')
    tprime <- 3650

    pop_child <- nrow(carriers_t %>% filter(agegroup == 0))
    pop_adult <- nrow(carriers_t %>% filter(agegroup == 1))
    pop_elder <- nrow(carriers_t %>% filter(agegroup == 2))


    ### create disinc df

    country_inv <- read.csv(paste("Q:/Technical/Python/", countrynam, "localinv.csv", sep = ""))

    disinc_t <- full_join(country_inv, output_t)

    disinc_t <- disinc_t %>% mutate(disease_cas = carrprev*invasiveness*indivs*(tprime/365)) %>%
      # disease incidence and disease cases
      mutate(disease_inc = (disease_cas/indivs)*100000,
             disease_cas_lo = carrprev*invasiveness.low*indivs*(tprime/365),
             disease_cas_hi = carrprev*invasiveness.high*indivs*(tprime/365)) %>%
      mutate(disease_inc_lo = (disease_cas_lo/indivs)*100000,
             disease_inc_hi = (disease_cas_hi/indivs)*100000) #dis per 100,000 ppl

    # overall carriage prevalence and IPD by age group
    carrprevchildsero3 <- apply(carriers_t[carriers_t$agegroup == 0, 2:70], 2, sum)/pop_child
    carrprevadultsero3 <- apply(carriers_t[carriers_t$agegroup == 1, 2:70], 2, sum)/pop_adult
    carrpreveldersero3 <- apply(carriers_t[carriers_t$agegroup == 2, 2:70], 2, sum)/pop_elder

    disinc_t$carrprevchild <- carrprevchildsero3
    disinc_t$carrprevadult <- carrprevadultsero3
    disinc_t$carrprevelder <- carrpreveldersero3

    disinc_t$X <- as.integer(runif(nrow(disinc_t), 0,100)) # random numbers so we can easily join w other df

    disinc <- full_join(disinc, disinc_t)

    ## create agegrpstats df

    # disease cases by serotype
    disease_child <- disinc_t$carrprevchild*disinc_t$invasiveness*(tprime/365)*pop_child
    disease_adult <- disinc_t$carrprevadult*disinc_t$invasiveness*(tprime/365)*pop_adult
    disease_elder <- disinc_t$carrprevelder*disinc_t$invasiveness*(tprime/365)*pop_elder

    disease_child_lo <- disinc_t$carrprevchild*disinc_t$invasiveness.low*(tprime/365)*pop_child
    disease_adult_lo <- disinc_t$carrprevadult*disinc_t$invasiveness.low*(tprime/365)*pop_adult
    disease_elder_lo <- disinc_t$carrprevelder*disinc_t$invasiveness.low*(tprime/365)*pop_elder
    disease_child_hi <- disinc_t$carrprevchild*disinc_t$invasiveness.high*(tprime/365)*pop_child
    disease_adult_hi <- disinc_t$carrprevadult*disinc_t$invasiveness.high*(tprime/365)*pop_adult
    disease_elder_hi <- disinc_t$carrprevelder*disinc_t$invasiveness.high*(tprime/365)*pop_elder

    # overall IPD incidence by age group
    overalldis <- (sum(disease_child, disease_adult, disease_elder)/(indivs*(tprime/365)))*100000
    overalldischild <- (sum(disease_child)/(pop_child*(tprime/365)))*100000
    overalldisadult <- (sum(disease_adult)/(pop_adult*(tprime/365)))*100000
    overalldiselder <- (sum(disease_elder)/(pop_elder*(tprime/365)))*100000

    overalldis_lo <- (sum(disease_child_lo, disease_adult_lo, disease_elder_lo)/(indivs*(tprime/365)))*100000
    overalldischild_lo <- (sum(disease_child_lo)/(pop_child*(tprime/365)))*100000
    overalldisadult_lo <- (sum(disease_adult_lo)/(pop_adult*(tprime/365)))*100000
    overalldiselder_lo <- (sum(disease_elder_lo)/(pop_elder*(tprime/365)))*100000
    overalldis_hi <- (sum(disease_child_hi + disease_adult_hi + disease_elder_hi)/(indivs*(tprime/365)))*100000
    overalldischild_hi <- (sum(disease_child_hi)/(pop_child*(tprime/365)))*100000
    overalldisadult_hi <- (sum(disease_adult_hi)/(pop_adult*(tprime/365)))*100000
    overalldiselder_hi <- (sum(disease_elder_hi)/(pop_elder*(tprime/365)))*100000

    overallcarrprev <- sum(carriers_t$susc < 2)/indivs
    overallcarrprevchild <- sum(carriers_t[carriers_t$agegroup == 0,]$susc < 2)/pop_child
    overallcarrprevadult <- sum(carriers_t[carriers_t$agegroup == 1,]$susc < 2)/pop_adult
    overallcarrprevelder <- sum(carriers_t[carriers_t$agegroup == 2,]$susc < 2)/pop_elder

    X <- c('carrprev', 'disease inc', 'disinc_lo', 'disinc_hi')
    overall <- c(overallcarrprev, overalldis, overalldis_lo, overalldis_hi)
    child <- c(overallcarrprevchild, overalldischild, overalldischild_lo, overalldischild_hi)
    adult <- c(overallcarrprevadult, overalldisadult, overalldisadult_lo, overalldisadult_hi)
    elder <- c(overallcarrprevelder, overalldiselder, overalldiselder_lo, overalldiselder_hi)

    agegrpstats_t <- data.frame(X, overall, child, adult, elder)
    agegrpstats_t$simnum <- i
    agegrpstats_t$vaxperiod <- 'postPCV30'
    agegrpstats_t$country <- countrynam
    agegrpstats_t$description <- unique(output_t$description)
    agegrpstats_t$vaxstrat <- unique(output_t$vaxstrat)

    agegrpstats <- full_join(agegrpstats, agegrpstats_t)
  }

  # re-calculate IRR for adults because it was written wrong in the simulations
  for (i in unique(IRR$simnum)) {
    agegrpstats_t <- agegrpstats %>% filter(simnum == i)
    overalldisadult <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'prevacc' &
                                             agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult2 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV13' &
                                             agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult3 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV20' &
                                              agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult4 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV30' &
                                              agegrpstats_t$X == 'disease inc'),]$adult
    post13_adult <- overalldisadult2/overalldisadult
    post20_adult <- overalldisadult3/overalldisadult2
    post30_adult <- overalldisadult4/overalldisadult3
    all_adult <- overalldisadult4/overalldisadult
    IRR[which(IRR$simnum == i & IRR$X == 'post13 vs pre' & IRR$agegrp == 'adults'),]$IRR <- post13_adult
    IRR[which(IRR$simnum == i & IRR$X == 'post20 vs post13' & IRR$agegrp == 'adults'),]$IRR <- post20_adult
    IRR[which(IRR$simnum == i & IRR$X == 'post30 vs post20' & IRR$agegrp == 'adults'),]$IRR <- post30_adult
    IRR[which(IRR$simnum == i & IRR$X == 'post30 vs pre' & IRR$agegrp == 'adults'),]$IRR <- all_adult
  }

  #re-estimate IRR for those with zero cases of disease that threw division by 0 error by adding 10e-8 cases
  sims_zeroIRR <- setdiff(agegrpstats$simnum, IRR$simnum)
  sims_disinc <- unique((disinc %>% filter(vaxperiod == 'postPCV30'))$simnum)
  sims_zeroIRR <- sims_zeroIRR[which(sims_zeroIRR %in% sims_disinc)]

  for (i in sims_zeroIRR) {
    # print(i)
    agegrpstats_t <- agegrpstats %>% filter(simnum == i)
    disinc_t <- disinc %>% filter(simnum == i)
    carriers_t <- carrierstmax %>% filter(simnum == i)

    # prevax
    disinc_t1 <- disinc_t %>% filter(vaxperiod == 'prevacc')
    carriers_t1 <- carriers_t %>% filter(vaxperiod == 'prevacc')
    pop_child <- nrow(carriers_t1 %>% filter(agegroup == 0))
    pop_adult <- nrow(carriers_t1 %>% filter(agegroup == 1))
    pop_elder <- nrow(carriers_t1 %>% filter(agegroup == 2))
    stop_t <- 7300

    # disease_child <- disinc_t1['carrprevchild']*disinc_t1['invasiveness']*(stop_t/365)*pop_child
    # disease_adult <- disinc_t1['carrprevadult']*disinc_t1['invasiveness']*(stop_t/365)*pop_adult
    # disease_elder <- disinc_t1['carrprevelder']*disinc_t1['invasiveness']*(stop_t/365)*pop_elder
    
    disease_child <- disinc_t1$carrprevchild*disinc_t1$invasiveness*(stop_t/365)*pop_child
    disease_adult <- disinc_t1$carrprevadult*disinc_t1$invasiveness*(stop_t/365)*pop_adult
    disease_elder <- disinc_t1$carrprevelder*disinc_t1$invasiveness*(stop_t/365)*pop_elder

    disease_all <- c(sum(disease_child), sum(disease_adult), sum(disease_elder))

    # postPCV13
    disinc_t2 <- disinc_t %>% filter(vaxperiod == 'postPCV13')
    carriers_t2 <- carriers_t %>% filter(vaxperiod == 'postPCV13')
    pop_child2 <- nrow(carriers_t2 %>% filter(agegroup == 0))
    pop_adult2 <- nrow(carriers_t2 %>% filter(agegroup == 1))
    pop_elder2 <- nrow(carriers_t2 %>% filter(agegroup == 2))
    stop_t <- 3650

    # disease_child2 <- disinc_t2['carrprevchild']*disinc_t2['invasiveness']*(stop_t/365)*pop_child2
    # disease_adult2 <- disinc_t2['carrprevadult']*disinc_t2['invasiveness']*(stop_t/365)*pop_adult2
    # disease_elder2 <- disinc_t2['carrprevelder']*disinc_t2['invasiveness']*(stop_t/365)*pop_elder2
    
    disease_child2 <- disinc_t2$carrprevchild*disinc_t2$invasiveness*(stop_t/365)*pop_child2
    disease_adult2 <- disinc_t2$carrprevadult*disinc_t2$invasiveness*(stop_t/365)*pop_adult2
    disease_elder2 <- disinc_t2$carrprevelder*disinc_t2$invasiveness*(stop_t/365)*pop_elder2

    disease_all2 <- c(sum(disease_child2), sum(disease_adult2), sum(disease_elder2))

    # postPCV20
    disinc_t3 <- disinc_t %>% filter(vaxperiod == 'postPCV20')
    carriers_t3 <- carriers_t %>% filter(vaxperiod == 'postPCV20')
    pop_child3 <- nrow(carriers_t3 %>% filter(agegroup == 0))
    pop_adult3 <- nrow(carriers_t3 %>% filter(agegroup == 1))
    pop_elder3 <- nrow(carriers_t3 %>% filter(agegroup == 2))

    # disease_child3 <- disinc_t3['carrprevchild']*disinc_t3['invasiveness']*(stop_t/365)*pop_child3
    # disease_adult3 <- disinc_t3['carrprevadult']*disinc_t3['invasiveness']*(stop_t/365)*pop_adult3
    # disease_elder3 <- disinc_t3['carrprevelder']*disinc_t3['invasiveness']*(stop_t/365)*pop_elder3
    
    disease_child3 <- disinc_t3$carrprevchild*disinc_t3$invasiveness*(stop_t/365)*pop_child3
    disease_adult3 <- disinc_t3$carrprevadult*disinc_t3$invasiveness*(stop_t/365)*pop_adult3
    disease_elder3 <- disinc_t3$carrprevelder*disinc_t3$invasiveness*(stop_t/365)*pop_elder3

    disease_all3 <- c(sum(disease_child3), sum(disease_adult3), sum(disease_elder3))

    # postPCV30
    disinc_t4 <- disinc_t %>% filter(vaxperiod == 'postPCV30')
    carriers_t4 <- carriers_t %>% filter(vaxperiod == 'postPCV30')
    pop_child4 <- nrow(carriers_t4 %>% filter(agegroup == 0))
    pop_adult4 <- nrow(carriers_t4 %>% filter(agegroup == 1))
    pop_elder4 <- nrow(carriers_t4 %>% filter(agegroup == 2))

    # disease_child4 <- disinc_t4['carrprevchild']*disinc_t4['invasiveness']*(stop_t/365)*pop_child4
    # disease_adult4 <- disinc_t4['carrprevadult']*disinc_t4['invasiveness']*(stop_t/365)*pop_adult4
    # disease_elder4 <- disinc_t4['carrprevelder']*disinc_t4['invasiveness']*(stop_t/365)*pop_elder4
    
    disease_child4 <- disinc_t4$carrprevchild*disinc_t4$invasiveness*(stop_t/365)*pop_child4
    disease_adult4 <- disinc_t4$carrprevadult*disinc_t4$invasiveness*(stop_t/365)*pop_adult4
    disease_elder4 <- disinc_t4$carrprevelder*disinc_t4$invasiveness*(stop_t/365)*pop_elder4

    disease_all4 <- c(sum(disease_child4), sum(disease_adult4), sum(disease_elder4))

    sumdisdf <- rbind(disease_all, disease_all2, disease_all3, disease_all4)
    colnames(sumdisdf) <- c('children', 'adults', 'elder')
    rownames(sumdisdf) <- c('prevax', 'postPCV13', 'postPCV20', 'postPCV30')

    sumdisdf[sumdisdf == 0] <- 1e-8

    # children
    overalldischild <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'prevacc' &
                                             agegrpstats_t$X == 'disease inc'),]$child
    overalldischild2 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV13' &
                                              agegrpstats_t$X == 'disease inc'),]$child
    overalldischild3 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV20' &
                                              agegrpstats_t$X == 'disease inc'),]$child
    overalldischild4 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV30' &
                                              agegrpstats_t$X == 'disease inc'),]$child
    post13_child <- overalldischild2/overalldischild
    post20_child <- overalldischild3/overalldischild2
    post30_child <- overalldischild4/overalldischild3
    all_child <- overalldischild4/overalldischild
    IRR_child <- c(post13_child, post20_child, post30_child, all_child)

    # stdev
    stdev_post13_child <- sqrt((1/sumdisdf['prevax', 'children']) + (1/sumdisdf['postPCV13', 'children']))
    stdev_post20_child <- sqrt((1/sumdisdf['postPCV13', 'children']) + (1/sumdisdf['postPCV20', 'children']))
    stdev_post30_child <- sqrt((1/sumdisdf['postPCV20', 'children']) + (1/sumdisdf['postPCV30', 'children']))
    stdev_all_child <- sqrt((1/sumdisdf['prevax', 'children']) + (1/sumdisdf['postPCV30', 'children']))

    # bounds 95% CI
    post13_child_bounds <- c(exp(log(post13_child)-(1.96*stdev_post13_child)),
                                    exp(log(post13_child)+(1.96*stdev_post13_child)))
    post20_child_bounds <- c(exp(log(post20_child)-(1.96*stdev_post20_child)),
                                    exp(log(post20_child)+(1.96*stdev_post20_child)))
    post30_child_bounds <- c(exp(log(post30_child)-(1.96*stdev_post30_child)),
                                    exp(log(post30_child)+(1.96*stdev_post30_child)))
    all_child_bounds <- c(exp(log(all_child)-(1.96*stdev_all_child)),
                                 exp(log(all_child)+(1.96*stdev_all_child)))

    IRR_bounds_child = rbind(post13_child_bounds, post20_child_bounds, post30_child_bounds, all_child_bounds)
    colnames(IRR_bounds_child) <- c('low', 'high')
    rownames(IRR_bounds_child) <- c('post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre')
    IRR_children <- data.frame(IRR_bounds_child)
    IRR_children$IRR <- IRR_child
    IRR_children$agegrp <- 'children'
    IRR_children$simnum <- i
    IRR_children$vaxperiod <- 'postPCV30' # irrelevant
    IRR_children$country <- 'fra'
    IRR_children$description <- simparams[which(simparams$simnum == simnum[i]),]$description
    IRR_children$vaxstrat <- 'local_inv'
    IRR_children <- tibble::rownames_to_column(IRR_children, var = 'X')

    #adults
    overalldisadult <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'prevacc' &
                                             agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult2 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV13' &
                                              agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult3 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV20' &
                                              agegrpstats_t$X == 'disease inc'),]$adult
    overalldisadult4 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV30' &
                                              agegrpstats_t$X == 'disease inc'),]$adult
    post13_adult <- overalldisadult2/overalldisadult
    post20_adult <- overalldisadult3/overalldisadult2
    post30_adult <- overalldisadult4/overalldisadult3
    all_adult <- overalldisadult4/overalldisadult
    IRR_adult <- c(post13_adult, post20_adult, post30_adult, all_adult)

    # stdev
    stdev_post13_adult <- sqrt((1/sumdisdf['prevax', 'adults']) + (1/sumdisdf['postPCV13', 'adults']))
    stdev_post20_adult <- sqrt((1/sumdisdf['postPCV13', 'adults']) + (1/sumdisdf['postPCV20', 'adults']))
    stdev_post30_adult <- sqrt((1/sumdisdf['postPCV20', 'adults']) + (1/sumdisdf['postPCV30', 'adults']))
    stdev_all_adult <- sqrt((1/sumdisdf['prevax', 'adults']) + (1/sumdisdf['postPCV30', 'adults']))

    # bounds 95% CI
    post13_adult_bounds <- c(exp(log(post13_adult)-(1.96*stdev_post13_adult)),
                             exp(log(post13_adult)+(1.96*stdev_post13_adult)))
    post20_adult_bounds <- c(exp(log(post20_adult)-(1.96*stdev_post20_adult)),
                             exp(log(post20_adult)+(1.96*stdev_post20_adult)))
    post30_adult_bounds <- c(exp(log(post30_adult)-(1.96*stdev_post30_adult)),
                             exp(log(post30_adult)+(1.96*stdev_post30_adult)))
    all_adult_bounds <- c(exp(log(all_adult)-(1.96*stdev_all_adult)),
                          exp(log(all_adult)+(1.96*stdev_all_adult)))

    IRR_bounds_adult = rbind(post13_adult_bounds, post20_adult_bounds, post30_adult_bounds, all_adult_bounds)
    colnames(IRR_bounds_adult) <- c('low', 'high')
    rownames(IRR_bounds_adult) <- c('post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre')
    IRR_adults <- data.frame(IRR_bounds_adult)
    IRR_adults$IRR <- IRR_adult
    IRR_adults$agegrp <- 'adults'
    IRR_adults$simnum <- i
    IRR_adults$vaxperiod <- 'postPCV30' # irrelevant
    IRR_adults$country <- 'fra'
    IRR_adults$description <- simparams[which(simparams$simnum == simnum[i]),]$description
    IRR_adults$vaxstrat <- 'local_inv'
    IRR_adults <- tibble::rownames_to_column(IRR_adults, var = 'X')

    # elder
    overalldiselder <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'prevacc' &
                                             agegrpstats_t$X == 'disease inc'),]$elder
    overalldiselder2 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV13' &
                                              agegrpstats_t$X == 'disease inc'),]$elder
    overalldiselder3 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV20' &
                                              agegrpstats_t$X == 'disease inc'),]$elder
    overalldiselder4 <- agegrpstats_t[which(agegrpstats_t$vaxperiod == 'postPCV30' &
                                              agegrpstats_t$X == 'disease inc'),]$elder
    overalldiselder_vec <- c(overalldiselder, overalldiselder2, overalldiselder3, overalldiselder4)
    overalldiselder_vec[overalldiselder_vec == 0] <- 1e-8
    post13_elder <- overalldiselder_vec[2]/overalldiselder_vec[1]
    post20_elder <- overalldiselder_vec[3]/overalldiselder_vec[2]
    post30_elder <- overalldiselder_vec[4]/overalldiselder_vec[3]
    all_elder <- overalldiselder_vec[4]/overalldiselder_vec[1]
    IRR_elder <- c(post13_elder, post20_elder, post30_elder, all_elder)

    # stdev
    stdev_post13_elder <- sqrt((1/sumdisdf['prevax', 'elder']) + (1/sumdisdf['postPCV13', 'elder']))
    stdev_post20_elder <- sqrt((1/sumdisdf['postPCV13', 'elder']) + (1/sumdisdf['postPCV20', 'elder']))
    stdev_post30_elder <- sqrt((1/sumdisdf['postPCV20', 'elder']) + (1/sumdisdf['postPCV30', 'elder']))
    stdev_all_elder <- sqrt((1/sumdisdf['prevax', 'elder']) + (1/sumdisdf['postPCV30', 'elder']))

    # bounds 95% CI
    post13_elder_bounds <- c(exp(log(post13_elder)-(1.96*stdev_post13_elder)),
                             exp(log(post13_elder)+(1.96*stdev_post13_elder)))
    post20_elder_bounds <- c(exp(log(post20_elder)-(1.96*stdev_post20_elder)),
                             exp(log(post20_elder)+(1.96*stdev_post20_elder)))
    post30_elder_bounds <- c(exp(log(post30_elder)-(1.96*stdev_post30_elder)),
                             exp(log(post30_elder)+(1.96*stdev_post30_elder)))
    all_elder_bounds <- c(exp(log(all_elder)-(1.96*stdev_all_elder)),
                          exp(log(all_elder)+(1.96*stdev_all_elder)))

    IRR_bounds_elder = rbind(post13_elder_bounds, post20_elder_bounds, post30_elder_bounds, all_elder_bounds)
    colnames(IRR_bounds_elder) <- c('low', 'high')
    rownames(IRR_bounds_elder) <- c('post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre')
    IRR_elders <- data.frame(IRR_bounds_elder)
    IRR_elders$IRR <- IRR_elder
    IRR_elders$agegrp <- 'elders'
    IRR_elders$simnum <- i
    IRR_elders$vaxperiod <- 'postPCV30' # irrelevant
    IRR_elders$country <- 'fra'
    IRR_elders$description <- simparams[which(simparams$simnum == simnum[i]),]$description
    IRR_elders$vaxstrat <- 'local_inv'
    IRR_elders <- tibble::rownames_to_column(IRR_elders, var = 'X')

    IRR_full <- bind_rows(IRR_children, IRR_adults, IRR_elders)
    IRR <- bind_rows(IRR, IRR_full)

    #print(i)

  }
  
  if (vaxstrat == 'local_inv') {
    IRR_locinv <<- IRR
    agegrpstats_locinv <<- agegrpstats
    disinc_locinv <<- disinc
    carrierstmax_locinv <<- carrierstmax
  }

  IRR_agegrps <- ggplot(IRR) +
    geom_errorbar(aes(x = X, ymin = low, ymax = high, colour = agegrp, group = agegrp),
                  width = 0.1, position = position_dodge(width = 0.3)) +
    geom_point(aes(x = X, y = IRR, colour = agegrp, group = agegrp),
               position = position_dodge(width = 0.3)) +
    #scale_y_continuous(trans = 'log10') +
    scale_color_viridis(discrete = T)

  IRR_description_ch <- ggplot(IRR %>% filter(agegrp == 'children')) +
    geom_errorbar(aes(x = X, ymin = low, ymax = high, colour = description, group = description),
                  width = 0.1, position = position_dodge(width = 0.3)) +
    geom_point(aes(x = X, y = IRR, colour = description, group = description),
               position = position_dodge(width = 0.3)) +
    scale_y_continuous(trans = 'log10') +
    scale_color_viridis(discrete = T)

  meta_df <- list()
  dflen <- 1

  for (i in unique(IRR$description)[grepl('5k', unique(IRR$description))]) { # only 5k simulations
    for (j in unique(IRR$X)) { # vax periods
      for (k in unique(IRR$agegrp)) { # each age group
        print(paste(i, j, k, sep = " - "))
        r <- pooled_est(df = IRR, des = i, vaxperiod2 = j, agegroup = k)
        meta_df[[dflen]] <- c(i, j, k, r$TE.fixed, r$lower.fixed, r$upper.fixed, r$TE.random, r$lower.random, r$upper.random)
        dflen <- dflen + 1
      }
    }
  }

  meta_df2 <- data.frame(do.call("rbind",meta_df))
  colnames(meta_df2) <- c('description', 'vaxperiod', 'agegrp',
                          'fixed', 'fixed.lo', 'fixed.hi',
                          'random', 'random.lo', 'random.hi')
  meta_df2$fixed <- as.character(meta_df2$fixed)
  meta_df2$fixed.lo <- as.character(meta_df2$fixed.lo)
  meta_df2$fixed.hi <- as.character(meta_df2$fixed.hi)
  meta_df2$random <- as.character(meta_df2$random)
  meta_df2$random.lo <- as.character(meta_df2$random.lo)
  meta_df2$random.hi <- as.character(meta_df2$random.hi)

  meta_df2$fixed <- as.numeric(meta_df2$fixed)
  meta_df2$fixed.lo <- as.numeric(meta_df2$fixed.lo)
  meta_df2$fixed.hi <- as.numeric(meta_df2$fixed.hi)
  meta_df2$random <- as.numeric(meta_df2$random)
  meta_df2$random.lo <- as.numeric(meta_df2$random.lo)
  meta_df2$random.hi <- as.numeric(meta_df2$random.hi)

  meta_df2$agegrp <- factor(meta_df2$agegrp, levels = c('all', 'children', 'adults', 'elders'),
                            labels = c('All', 'Children', 'Adults', 'Elderly'))
  meta_df2$vaxperiod <- factor(meta_df2$vaxperiod, levels = c('post13 vs pre', 'post20 vs post13',
                                                              'post30 vs post20', 'post30 vs pre'),
                            labels = c('PCV13', 'PCV20', 'PCV30', 'All PCVs'))
  meta_df2$description <- factor(meta_df2$description,
                                 levels = c("default5k",
                                            "eps5k_lo", "eps5k_hi",
                                            "sigma5k_lo", "sigma5k_hi",
                                            "comp5k_lo", "comp5k_hi",
                                            "vaxeff5k_lo", "vaxeff5k_hi",
                                            "vaxdur5k_lo", "vaxdur5k_hi"),
                                 labels = c('Default', 'Epsilon = 0.25', 'Epsilon = 0.5',
                                            'Sigma = 0.6', 'Sigma = 0.8',
                                            'Competition = 0.2', 'Competition = 0.5',
                                            'Vaccine efficacy = 0.75', 'Vaccine efficacy = 0.85',
                                            'Vaccine duration = 5 years', 'Vaccine duration = 25 years'))
  meta_df2$vaxstrategy <- vaxstrat
  meta_df2$country <- countrynam

  # r <- list()
  # q <- 1
  # for (i in unique(meta_df2$description)) {
  #   s <- pl(i, metadataframe = meta_df2)
  #   r[[q]] <- s
  #   q <- q+1
  # }
  # 
  # IRR_plot <- ggarrange(r[[1]], r[[11]], r[[2]], "", r[[9]], r[[10]], "", r[[5]],
  #           r[[6]], "", r[[7]], r[[8]], "", r[[3]], r[[4]],
  #           ncol = 3, nrow = 5,
  #           common.legend = T, legend = 'bottom')
  #ggsave(IRR_plot, file = 'fra_inv_allIRR.pdf', height = 10, width = 9)
  return(meta_df2)
}

## fra
fra_inv_meta <- IRRfun(fra_inv_file)
fra_IRR_locinv <- IRR_locinv
fra_agegrpstats_locinv <- agegrpstats_locinv
fra_disinc_locinv <- disinc_locinv
fra_carrierstmax_locinv <- carrierstmax_locinv

fra_dis_meta <- IRRfun(fra_dis_file)
fra_dis_agegrpstats <- agegrpstats
fra_dis_VTsero <- VTsero
fra_glob_meta <- IRRfun(fra_glob_file)
fra_glob_agegrpstats <- agegrpstats

## usa
usa_inv_meta <- IRRfun(usa_inv_file)
usa_IRR_locinv <- IRR_locinv
usa_agegrpstats_locinv <- agegrpstats_locinv
usa_disinc_locinv <- disinc_locinv
usa_carrierstmax_locinv <- carrierstmax_locinv

usa_dis_meta <- IRRfun(usa_dis_file)
usa_dis_agegrpstats <- agegrpstats
usa_glob_meta <- IRRfun(usa_glob_file)
usa_glob_agegrpstats <- agegrpstats

####################################################################################################
### Carriage prevalence analysis
####################################################################################################
# fra
fra_agegrpstats <- bind_rows(fra_agegrpstats_locinv, fra_dis_agegrpstats, fra_glob_agegrpstats)
fra_carrprev5k <- fra_agegrpstats[grepl('5k', fra_agegrpstats$description),] %>% filter(X == 'carrprev') # only 5k sims
fra_carrprev5k$vaxstrat <- factor(fra_carrprev5k$vaxstrat, 
                                  levels = c('dis_inc', 'local_inv', 'glob_inv'),
                                  labels = c('Disease incidence', 'Local invasiveness', 'Global invasiveness'))
counts_fr <- fra_carrprev5k %>% group_by(description, vaxstrat, vaxperiod) %>% summarise(count = n())


# usa

usa_agegrpstats <- bind_rows(usa_agegrpstats_locinv, usa_dis_agegrpstats, usa_glob_agegrpstats)
usa_carrprev5k <- usa_agegrpstats[grepl('5k', usa_agegrpstats$description),] %>% filter(X == 'carrprev') # only 5k sims
usa_carrprev5k$vaxstrat <- factor(usa_carrprev5k$vaxstrat, 
                                  levels = c('dis_inc', 'local_inv', 'glob_inv'),
                                  labels = c('Disease incidence', 'Local invasiveness', 'Global invasiveness'))
counts_us <- usa_carrprev5k %>% group_by(description, vaxstrat, vaxperiod) %>% summarise(count = n())

lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}
upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

temp <- fra_carrprev5k %>% 
  pivot_longer(cols = c('overall', 'child', 'adult', 'elder'), 
               names_to = 'agegrp', 
               values_to = 'carrprev') %>%
  distinct()
temp2 <- usa_carrprev5k %>%
  pivot_longer(cols = c('overall', 'child', 'adult', 'elder'), 
               names_to = 'agegrp', 
               values_to = 'carrprev') %>%
  distinct()

fra_carrprev_avg <- temp %>% group_by(description, vaxperiod, vaxstrat, agegrp) %>% 
  summarise(smean = mean(carrprev, na.rm = TRUE),
            ssd = sd(carrprev, na.rm = TRUE),
            count = n()) %>%
  mutate(se = ssd / sqrt(count),
         lower_ci = lower_ci(smean, se, count),
         upper_ci = upper_ci(smean, se, count),
         country = 'fra')


usa_carrprev_avg <- temp2 %>% group_by(description, vaxperiod, vaxstrat, agegrp) %>% 
  summarise(smean = mean(carrprev, na.rm = TRUE),
            ssd = sd(carrprev, na.rm = TRUE),
            count = n()) %>%
  mutate(se = ssd / sqrt(count),
         lower_ci = lower_ci(smean, se, count),
         upper_ci = upper_ci(smean, se, count),
         country = 'usa')

carrprev_avg <- bind_rows(fra_carrprev_avg, usa_carrprev_avg)

carrprev_avg$description <- factor(carrprev_avg$description,
                               levels = c("default5k",
                                          "eps5k_lo", "eps5k_hi",
                                          "sigma5k_lo", "sigma5k_hi",
                                          "comp5k_lo", "comp5k_hi",
                                          "vaxeff5k_lo", "vaxeff5k_hi",
                                          "vaxdur5k_lo", "vaxdur5k_hi"),
                               labels = c('Default', 'Epsilon = 0.25', 'Epsilon = 0.5',
                                          'Sigma = 0.6', 'Sigma = 0.8',
                                          'Competition = 0.2', 'Competition = 0.5',
                                          'Vaccine efficacy = 0.75', 'Vaccine efficacy = 0.85',
                                          'Vaccine duration = 5 years', 'Vaccine duration = 25 years'))
carrprev_avg$agegrp <- factor(carrprev_avg$agegrp,
                             levels = c('overall', 'child', 'adult', 'elder'),
                             labels = c('Overall', 'Children', 'Adults', 'Elderly'))

### table with default scenario carr prev averages 
newtemp <- carrprev_avg %>% filter(description == 'Default')
tab6 <- newtemp %>% 
  mutate(vals = paste(round(smean*100, digits = 1), '% (', 
                      round(100*lower_ci, digits = 1), '% - ', 
                      round(100*upper_ci, digits = 1), '%)', sep = '')) %>%
  select(vaxperiod, vaxstrat, agegrp, vals, country) %>%
  pivot_wider(names_from = c(agegrp, vaxperiod), names_sep = '_', values_from = vals) 
tab6 <- tab6[,c(1,2,3,7,5,4,6,11,9,8,10,15,13,12,14)]
setwd("Q:/Technical/Python/figs for paper")
write.csv(tab6, 'tab6_carrprev_postPCV_allstratallage.csv')

# prevax
carrprevplot <- function(agegroup) {
  carrprev_avg$country <- factor(carrprev_avg$country, levels = c('fra', 'usa'), labels = c('France', 'USA'))
  ggplot(carrprev_avg %>% filter(agegrp == agegroup, vaxperiod == 'prevacc'), 
         group = country) + 
  geom_errorbar(aes(x = description, ymin = lower_ci*100, ymax = upper_ci*100, colour = country), 
                width = 0.1, position = position_dodge(0.2)) + 
  geom_point(aes(x = description, y = smean*100, colour = country), position = position_dodge(0.2)) + 
  labs(x = "", y = "Carriage prevalence (%)", colour = '') +
  theme_bw() + theme_minimal() + ggtitle(agegroup) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis(discrete = T)
}

# postvax
carrprevplot_post <- function(agegroup, vaccinperiod, countrrry) {
  ggplot(carrprev_avg %>% filter(vaxperiod == vaccinperiod, agegrp == agegroup, country == countrrry), 
         aes(group = vaxstrat)) + 
    geom_errorbar(aes(x = description, 
                      ymin = lower_ci*100, 
                      ymax = upper_ci*100, 
                      colour = vaxstrat), 
                  width = 0.1, position = position_dodge(0.3)) + 
    geom_point(aes(x = description, y = smean*100, colour = vaxstrat), position = position_dodge(0.3)) + 
    labs(x = "", y = 'Carriage prevalence (%)', colour = "") + 
    theme_bw() + theme_minimal() + ggtitle(agegroup) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_color_viridis(discrete = T)
}

# pre-PCV
carrprev_pre_plots <- lapply(unique(carrprev_avg$agegrp), function(x) carrprevplot(agegroup = x))
ggarrange(carrprev_pre_plots[[4]], carrprev_pre_plots[[2]], carrprev_pre_plots[[1]], carrprev_pre_plots[[3]],
          nrow = 2, ncol = 2, common.legend = T, legend = 'bottom') # pdf 6 x 8

# post-PCV13
fra_carrprev_post13_plots <- lapply(unique(carrprev_avg$agegrp), 
                                function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV13', countrrry = 'fra'))
fra_post13_carrprev_pl <- ggarrange(fra_carrprev_post13_plots[[4]], fra_carrprev_post13_plots[[2]], 
                                    fra_carrprev_post13_plots[[1]], fra_carrprev_post13_plots[[3]],
          nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(fra_post13_carrprev_pl, top = text_grob("France - Post-hPCV13", face = "bold", size = 15)) # pdf 6x8

usa_carrprev_post13_plots <- lapply(unique(carrprev_avg$agegrp), 
                                    function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV13', countrrry = 'usa'))
usa_post13_carrprev_pl <- ggarrange(usa_carrprev_post13_plots[[4]], usa_carrprev_post13_plots[[2]], 
                                    usa_carrprev_post13_plots[[1]], usa_carrprev_post13_plots[[3]],
                                    nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(usa_post13_carrprev_pl, top = text_grob("USA - Post-hPCV13", face = "bold", size = 15)) # pdf 6x8

# post-PCV20
fra_carrprev_post20_plots <- lapply(unique(carrprev_avg$agegrp), 
                                    function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV20', countrrry = 'fra'))
fra_post20_carrprev_pl <- ggarrange(fra_carrprev_post20_plots[[4]], fra_carrprev_post20_plots[[2]], 
                                    fra_carrprev_post20_plots[[1]], fra_carrprev_post20_plots[[3]],
                                    nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(fra_post20_carrprev_pl, top = text_grob("France - Post-hPCV20", face = "bold", size = 15)) # pdf 6x8

usa_carrprev_post20_plots <- lapply(unique(carrprev_avg$agegrp), 
                                    function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV20', countrrry = 'usa'))
usa_post20_carrprev_pl <- ggarrange(usa_carrprev_post20_plots[[4]], usa_carrprev_post20_plots[[2]], 
                                    usa_carrprev_post20_plots[[1]], usa_carrprev_post20_plots[[3]],
                                    nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(usa_post20_carrprev_pl, top = text_grob("USA - Post-hPCV20", face = "bold", size = 15)) # pdf 6x8

# post-PCV30
fra_carrprev_post30_plots <- lapply(unique(carrprev_avg$agegrp), 
                                    function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV30', countrrry = 'fra'))
fra_post30_carrprev_pl <- ggarrange(fra_carrprev_post30_plots[[4]], fra_carrprev_post30_plots[[2]], 
                                    fra_carrprev_post30_plots[[1]], fra_carrprev_post30_plots[[3]],
                                    nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(fra_post30_carrprev_pl, top = text_grob("France - Post-hPCV30", face = "bold", size = 15)) # pdf 6x8

usa_carrprev_post30_plots <- lapply(unique(carrprev_avg$agegrp), 
                                    function(x) carrprevplot_post(agegroup = x, vaccinperiod = 'postPCV30', countrrry = 'usa'))
usa_post30_carrprev_pl <- ggarrange(usa_carrprev_post30_plots[[4]], usa_carrprev_post30_plots[[2]], 
                                    usa_carrprev_post30_plots[[1]], usa_carrprev_post30_plots[[3]],
                                    nrow = 2, ncol = 2, common.legend = T, legend = 'bottom')
annotate_figure(usa_post30_carrprev_pl, top = text_grob("USA - Post-hPCV30", face = "bold", size = 15)) # pdf 6x8


####################################################################################################
### Pooled IRR analysis
####################################################################################################
fra_meta <- bind_rows(fra_inv_meta, fra_dis_meta, fra_glob_meta)
fra_meta$vaxstrategy <- factor(fra_meta$vaxstrategy, levels = c('dis_inc', 'local_inv', 'glob_inv'),
                               labels = c('Disease incidence', 'Local invasiveness', 'Global invasiveness'))

### plot with vaccine vs pooled IRR by sensitivity analysis with shape as vax strategy and colour as age grp
r <- list()
q <- 1
for (i in unique(fra_meta$description)) {
  s <- pl(i, metadataframe = fra_meta)
  r[[q]] <- s
  q <- q+1
}

# immunity competition
IRR_plot <- ggarrange(r[[1]], r[[11]], r[[2]], "", r[[9]], r[[10]], "", r[[5]],
          r[[6]],
          ncol = 3, nrow = 3,
          common.legend = T, legend = 'bottom')
ggsave(IRR_plot, file = 'Q:/Technical/Python/figs for paper/fra_IRR_immcomp.pdf', height = 10, width = 9)

# vaccine
IRR_plot2 <- ggarrange(r[[1]], r[[7]], r[[8]], "", r[[3]], r[[4]],
                      ncol = 3, nrow = 2,
                      common.legend = T, legend = 'bottom')
ggsave(IRR_plot2, file = 'Q:/Technical/Python/figs for paper/fra_IRR_vax.pdf', height = 10, width = 9)

### plot with vaccine vs pooled IRR by age grp; colour is by sensitivity analysis 
r <- list()
q <- 1
for (i in unique(fra_meta$vaxstrategy)) {
  for (j in unique(fra_meta$agegrp)) {
    s <- pl2(j, metadataframe = fra_meta, vaxstrat = i)
    r[[q]] <- s
    q <- q+1
  }
  
}

plot_localinv <- ggarrange(r[[1]], r[[2]], r[[3]], r[[4]],
                      ncol = 2, nrow = 2,
                      common.legend = T, legend = 'right')
ggsave(plot_localinv, file = 'Q:/Technical/Python/figs for paper/fra_localinv_IRR.pdf', height = 10, width = 9)

plot_disinc <- ggarrange(r[[5]], r[[6]], r[[7]], r[[8]],
                      ncol = 2, nrow = 2,
                      common.legend = T, legend = 'right')
ggsave(plot_disinc, file = 'Q:/Technical/Python/figs for paper/fra_disinc_IRR.pdf', height = 10, width = 9)

plot_globalinv <- ggarrange(r[[9]], r[[10]], r[[11]], r[[12]],
                      ncol = 2, nrow = 2,
                      common.legend = T, legend = 'right')
ggsave(plot_globalinv, file = 'Q:/Technical/Python/popsims/fra_globalinv_IRR.pdf', height = 10, width = 9)

usa_inv_meta <- IRRfun(usa_inv_file)
usa_dis_meta <- IRRfun(usa_dis_file)
usa_glob_meta <- IRRfun(usa_glob_file)

####################################################################################################
### Which VT most frequent for each hPCV?
####################################################################################################

## for dis inc only: which VT most frequent for each PCV?
freqPCV13 <- count(VTsero, vt_1) %>% arrange(desc(n))
freqPCV13_VT <- freqPCV13$vt_1[1:13]
freqPCV20 <- count(VTsero, vt_2) %>% arrange(desc(n))
freqPCV20_VT <- freqPCV13$vt_1[1:20]
freqPCV30 <- count(VTsero, vt_3) %>% arrange(desc(n))
freqPCV30_VT <- freqPCV13$vt_1[1:30]

## how many simulations run of each model?
nu <- disinc %>% group_by(simnum, country, description, vaxstrat) %>% nest()
simnums <- table(nu$description)

## which sims didn't run?
full_fr <- c(1:223)
setdiff(full_fr, nu$simnum)

full_usa <- c(224:414)
setdiff(full_usa, nu$simnum) 