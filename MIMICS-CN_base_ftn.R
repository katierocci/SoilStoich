## Set working drive
rm(list = ls())
#setwd("/Users/wwieder/Will/git_repos_local/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)
library(deSolve)
library(tidyverse)
library(ggpattern)

#bring in RXEQ function
source("CN_RXEQ.R")
source("calc_Tpars.R")


###########################################
# MIMICS single point function - priming and desorption experiment - Figure 4
###########################################
experiment = c('Control','Priming','Priming+Acid')
for (x in 1:3)  {   # loop over exudation experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ

  Site = df$Site
  ANPP = df$ANPP
  TSOI = df$MAT
  CLAY = df$CLAY2 #5 or 55
  LIG = df$LIG
  CN = df$CN #45
  exud = 0.1

  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################


  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud,nUPmod=1) #>> Same example site input as used for stode


  #----------initialize pools---------------
  LIT_1  <<- 1e-4
  LIT_2  <<- 1e-4
  MIC_1  <<- 1e-4
  MIC_2  <<- 1e-4
  SOM_1  <<- 1e-4
  SOM_2  <<- 1e-4
  SOM_3  <<- 1e-4

  LIT_1_N  <<- 1e-4
  LIT_2_N  <<- 1e-4
  MIC_1_N  <<- 1e-4
  MIC_2_N  <<- 1e-4
  SOM_1_N  <<- 1e-4
  SOM_2_N  <<- 1e-4
  SOM_3_N  <<- 1e-4
  DIN      <<- 1e-4
  LeachingLoss      <<- 1e-4

  pools = c('LIT_m',  'LIT_s',  'MIC_r',  'MIC_K',  'SOM_p',  'SOM_c',  'SOM_a',
            'LIT_m_N','LIT_s_N','MIC_r_N','MIC_K_N','SOM_p_N','SOM_c_N','SOM_a_N',
            'DIN')

  # Update MIMICS pools
  #---------------------------------------------

  Ty    <<- c(LIT_1 = LIT_1, LIT_2 = LIT_2,
              MIC_1 = MIC_1, MIC_2 = MIC_2,
              SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3,
              LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N,
              MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N,
              SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
              DIN = DIN)

  test  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)

  # save results to dataframe
  if (x == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[x]
  } else if (x==2) {
    test2 = test
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[x]
    df_ss <- rbind(df_ss, df_exp)
  } else if (x==3) {
    test3 = test
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[x]
    df_ss <- rbind(df_ss, df_exp)
  }

}  # close experiment loop

test1[[1]][15]
test3[[1]][15]
test3
# MIC r:K
(test1[[1]][3]/test1[[1]][4])
(test2[[1]][3]/test2[[1]][4])

# bulk C:N
round(sum(test1[[1]][1:7])/sum(test1[[1]][8:14]),1)
round(sum(test2[[1]][1:7])/sum(test2[[1]][8:14]),1)
round(sum(test3[[1]][1:7])/sum(test3[[1]][8:14]),1)

#total N
sum(test1[[1]][8:15])
sum(test2[[1]][8:15])
sum(test3[[1]][8:15])

# MIC : TOTAL
sum(test1[[1]][3:4])/sum(test1[[1]][1:7])
sum(test2[[1]][3:4])/sum(test2[[1]][1:7])

# SOMP : SOM TOTAL
sum(test1[[1]][5])/sum(test1[[1]][5:7])
sum(test1[[1]][6])/sum(test1[[1]][5:7])

# MIC C:N
sum(test1[[1]][3:4])/sum(test1[[1]][10:11])
sum(test2[[1]][3:4])/sum(test2[[1]][10:11])

# Difference in N pools kg/m2)
(sum(test1[[1]][8:14]) - sum(test3[[1]][8:14])) * MICROtoECO * 1e3

#====================================
## Finished steady state calculations
## plot results
#====================================
Site = "Temperate deciduous forest"
df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)


df_ss %>%
  filter(as.numeric(pools) <= 7) %>%
  pivot_wider(id_cols = "pools", values_from = "y", names_from = "experiment") %>%
  mutate(Priming = Priming/Control,
         `Priming+Acid` = `Priming+Acid`/Control,
         Control = Control/Control) %>%
  pivot_longer(all_of(experiment), names_to = "Experiment") %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  filter(experiment != "Control") %>%
  ggplot(aes(x = pools, y = value, color = experiment))  +
  geom_point(size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
  scale_color_manual(values = c("lightskyblue", "darkgoldenrod2")) +
  geom_hline(yintercept = 1, color = "aquamarine4", linewidth = 2) +
  labs(y = "response Ratio (treatment/control)")+theme_bw()


df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  ggplot(aes(x=experiment, y=CNRatio,fill=experiment,show.legend = FALSE)) +
  geom_bar(stat = "identity",show.legend = FALSE) +
  scale_fill_manual(values = c("aquamarine4","lightskyblue", "darkgoldenrod2")) +
  coord_cartesian(ylim=c(8.5,10)) +
  labs(y = "Bulk C:N",
       title = paste("Bulk C:N: Exudation,",Site)) +theme_bw()


########################################
########################################
# call time series for elevated CO2 experiment - Figure 6a
# SS forward
#######################################
########################################

experiment = c('eCO2_baseline',
               'eCO2_prime','eCO2_prime_rootExud',
               'eCO2_prime_acid','eCO2_prime_acid_rootExud')
tx = c('baseline','prime','prime','prime+acid','prime+acid')
alloc = c('10%','10%','20%','10%','20%')
sim_year = 51 # Set number of days to sim forward
sim_days = 365*sim_year
exud = c(0.1,0.2)
eCO2_NPP = 1.2 #change to 1.0 to remove NPP effects
eCO2_CHEM = 1.1
nUPmod = 1

for (z in seq(1,length(experiment))) {
  # Grab steady state pools from stode output above
  if (z ==1 ) {
    ss = test1[[1]]
    x_adj = 1
  } else if (z <=3) {
    ss = test2[[1]]
    x_adj = 2
  } else {
    ss = test3[[1]]
    x_adj = 3
  }
  print(x_adj)
  # Create dataframe to store MIMICS pools over timesteps
  ss[["Nmin"]] <- 0
  ss[["Cover"]] <- 0
  MIMfwd = t(as.data.frame(ss))


  # Hourly model loop
  for(i in 2:(sim_days)){ # interval = hour

    # Recalc Tpars here, calls ftn in calc_Tpars.R
    #---------------------------------------------
    # e.g. 20% increase in NPP, starting in year 1
    exud_adj = exud[1]

    if (i < (1*365) ) {
      ANPP_adj = ANPP
      CN_adj = CN
    } else {
      ANPP_adj = ANPP * eCO2_NPP
      CN_adj = CN * eCO2_CHEM
      if  (z==3 || z==5)  {
        exud_adj = exud[2]
      }
    }
    #>> Same example site input as used for stode
    Tpars_mod = calc_Tpars(TSOI = TSOI, ANPP = ANPP_adj, CLAY = CLAY,
                           CN = CN_adj, LIG = LIG,x=x_adj,exud=exud_adj,
                           nUPmod=nUPmod)

    # Update MIMICS pools
    #---------------------------------------------
    step = CN_iter(t=NA, y=MIMfwd[i-1,], pars=Tpars_mod)
    pools_update = c(MIMfwd[i-1,1:15] + unlist(step)[1:15],
                     Nmin = unlist(step)[16],
                     Cover = unlist(step)[17])
    MIMfwd = rbind(MIMfwd, t(as.data.frame(pools_update)))
  }     # close daily (i) loop

  # Steady state forward hourly data
  if (z == 1){
    df_out = as.data.frame(MIMfwd)
    df_out$experiment <- rep(experiment[z],length(sim_days))
    df_out$tx    <- rep(tx[z],length(sim_days))
    df_out$alloc <- rep(alloc[z],length(sim_days))
    df_out$year  <- seq(0,sim_year,length.out=sim_days)
  } else {
    df_exp = as.data.frame(MIMfwd)
    df_exp$experiment <- rep(experiment[z],length(sim_days))
    df_exp$tx    <- rep(tx[z],length(sim_days))
    df_exp$alloc <- rep(alloc[z],length(sim_days))
    df_exp$year  <- seq(0,sim_year,length.out=sim_days)
    df_out <- rbind(df_out, df_exp)
  }

} # close z loop



df_out$bulkCN = rowSums(df_out[,1:7])/rowSums(df_out[,8:14])
df_out$TotalC = rowSums(df_out[,1:7]) * MICROtoECO
df_out$MicC = rowSums(df_out[,3:4]) * MICROtoECO
df_out$MicCN = rowSums(df_out[,3:4])/rowSums(df_out[,10:11])
df_out$MicC_O = df_out[,3]/df_out[,4]
df_out$Nmin  = df_out$Nmin * MICROtoECO
df_out$Cover  = df_out$Cover * MICROtoECO
df_out$experiment<- factor(df_out$experiment, levels = experiment)
df_out$alloc<- factor(df_out$alloc, levels = c('10%','20%'))


#plots of elevated CO2

d = df_out$year[2]
df_diff <- df_out %>%
  filter(year == sim_year | year == year[2]) %>%
  mutate(year = replace(year, year == d, "start")) %>%
  mutate(year = replace(year, year == 51.000000000, "end")) %>%
  pivot_wider(values_from = c(LIT_1,LIT_2, MIC_1, MIC_2, SOM_1, SOM_2, SOM_3, LIT_1_N,LIT_2_N, MIC_1_N, MIC_2_N, SOM_1_N, SOM_2_N, SOM_3_N, DIN, Nmin, bulkCN, TotalC, MicC, MicC_O, MicCN, Cover),
              names_from = year)
res <- df_diff[, grepl("_end", colnames(df_diff))]/ df_diff[, grepl("_start", colnames(df_diff))]
for ( col in 1:ncol(res)){
  colnames(res)[col] <-  sub("_end", "", colnames(res)[col])
}
df_front <- select(df_diff, experiment, tx, alloc)
df_diff <- cbind(df_front, res) #difference between 50 and 0 years in CO2
df_diff$tx <- factor(df_diff$tx, levels=c('baseline','prime', 'prime+acid'),
                     labels=c('Baseline', 'Priming', 'Priming+Desorption'))

#Figure 6a
level_order <- c('Bulk_CtoN', 'MicC', 'MicC_O', 'SOMc', 'N_min')
eCO2_pools.prim <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(SOMp = SOM_1) %>%
  mutate(SOMc = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, MicC, MicC_O, SOMc, N_min) %>%
  pivot_longer(cols= 4:8, values_to = "dif", names_to = "Pool") %>%
  filter(tx == "Priming") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) + ylim(0.8,1.3) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools", shape="Allocation to \nroot exudates") +
  theme_bw(base_size=16)
eCO2_pools.prim


#observed response ratios from Duke FACE experiment - Figure 6b
FACE_Obs <- read.csv("Drakeetal2011_DukeFACEdata.csv")
FACE_Obs <- FACE_Obs[1, ]
level_order <- c('Bulk_CtoN', 'MicC', 'Bac2Fun', 'fLF', 'N_min')
eCO2_obs <- FACE_Obs %>% select('Soil_CN_RR', 'MB_RR', 'BF_RR', 'fLF_RR', 'Nmin_RR') %>%
  mutate(Bulk_CtoN = as.numeric(Soil_CN_RR)) %>% mutate(MicC = as.numeric(MB_RR)) %>%
  mutate(Bac2Fun = as.numeric(BF_RR)) %>% mutate(fLF = as.numeric(fLF_RR)) %>% mutate(N_min = as.numeric(Nmin_RR)) %>%
  pivot_longer(cols = 6:10, names_to = 'Pool', values_to = 'RR') %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = RR))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(size=5, alpha=0.7, shape = 15) + scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3), limits = c(0.5, 1.3)) +
  labs(y = "Observed response ratio \nunder elevated CO2", x= " Pools") +
  theme_bw(base_size=16)
eCO2_obs


