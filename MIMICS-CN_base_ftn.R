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
# MIMICS single point function
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

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (kgC/m2)") +
  geom_text(x=1, y=1, label=paste("CN =",round(sum(test1[[1]][1:7])/sum(test1[[1]][8:14]),1))) +
  geom_text(x=2, y=1, label=paste("CN =",round(sum(test2[[1]][1:7])/sum(test2[[1]][8:14]),1))) +
  geom_text(x=3, y=1, label=paste("CN =",round(sum(test3[[1]][1:7])/sum(test3[[1]][8:14]),1)))


df_ss %>%
  filter(as.numeric(pools) <= 7) %>%
  pivot_wider(id_cols = "pools", values_from = "y", names_from = "experiment") %>%
  mutate(Priming = Priming/Control,
         `Priming+Acid` = `Priming+Acid`/Control,
         Control = Control/Control) %>%
  pivot_longer(all_of(experiment), names_to = "Experiment") %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  ggplot(aes(x = pools, y = value, color = experiment))  +
  geom_point(size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
  labs(y = "response ratio (tx/control)",
       title = paste("Exudation effects,",Site))


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
  coord_cartesian(ylim=c(6,12)) +
  labs(y = "Bulk C:N",
       title = paste("Bulk C:N: Exudation,",Site))


########################################
########################################
# call time series for elevated CO2 experiment...
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
eCO2_NPP = 1.0
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



# example pool change over simulation period
ggplot(df_out, aes(x=year,y = bulkCN,color=tx,linetype=alloc)) +
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotdash")) +
  xlab("Year")+
  labs(y = "Bulk C:N",
       title = paste("eCO2_NPP =",eCO2_NPP,", eCO2_chem =",eCO2_CHEM))

ggplot(df_out, aes(x=year,y = soilCN,color=tx,linetype=alloc)) +
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotdash")) +
  xlab("Year")+
  labs(y = "Soil C:N",
       title = paste("eCO2_NPP =",eCO2_NPP,", eCO2_chem =",eCO2_CHEM))


ggplot(df_out, aes(x=year,y = TotalC,color=tx,linetype=alloc)) +
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotdash")) +
  xlab("Year") +
  labs(y = "Total C",
     title = paste("eCO2_NPP =",eCO2_NPP,", eCO2_chem =",eCO2_CHEM))

ggplot(df_out, aes(x=year,y = MicC_O,color=tx,linetype=alloc)) +
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotdash")) +
  xlab("Year") +
  labs(y = "MIC C:O",
       title = paste("eCO2_NPP =",eCO2_NPP,", eCO2_chem =",eCO2_CHEM))

df_out %>%
  filter(year >= 1) %>%
  ggplot(aes(x=year,y = Nmin,color=tx,linetype=alloc)) +
  geom_line()+
  scale_linetype_manual(values = c("solid", "dotdash")) +
  xlab("Year")

ggplot(df_out, aes(x=year,y = (SOM_1/(SOM_1+SOM_2+SOM_3)),
                   color=experiment)) +
  geom_line()+
  ylab("MAOM fraction of total") +
  xlab("Year")

ggplot(df_out, aes(x=year,y = (SOM_1),
                   color=experiment)) +
  geom_line()+
  ylab("MAOM") +
  xlab("Year")

ggplot(df_out, aes(x=year,y = (SOM_2/(SOM_1+SOM_2+SOM_3)),
                   color=experiment)) +
  geom_line()+
  ylab("POM fraction of total") +
  xlab("Year")


ggplot(df_out, aes(x=year,y = (DIN),
                   color=experiment)) +
  geom_line()+
  xlab("Year")

ggplot(df_out, aes(x=year,y = (Cover),
                   color=experiment)) +
  geom_line()+
  xlab("Year")

df_out[1,]
nsteps = dim(df_out)[1]
df_out[nsteps,]

## Quick look at response ratios:
x <- df_out %>%
  filter(year==sim_year) #%>%
y <- df_out %>%
  filter(year==df_out$year[2]) #%>%
z = x[,1:16]/y[,1:16]
z
z2 = x[,22:26]/y[,22:26]
z2


# exudate effect kg/m2/d to g/m2/y to kg/ha
#exef_prime = 1e3*365*(4.177809e-05 - 4.168749e-05) * 1e4 *1e-3
#exef_acid = 1e3*365*(4.227151e-05 - 4.170081e-05)  * 1e4 *1e-3
#exef_acid
#exef_prime

# Nmin (kg/m2/d -> mgC/cm3/d (0-30 cm) -> ug/cm3/h)
4.167958e-05 / MICROtoECO * 1e3/24

#KTs plots down here
#Make a combined, just LQ, just NPP series of plots so we can see what is driving what?

#response ratio of end of simulation to second line of year (so not 0) - need to figure out the naming convention of the rows...

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
#plot of differences
eCO2_CN <- df_diff %>%
  mutate(Site == "Forest") %>%
  filter(tx == "Priming") %>%
  ggplot(aes(x = Site, y=bulkCN, shape =alloc))+ geom_point(size = 6, position = "jitter") +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +ylab("Change in bulk soil C:N after \n50 years of elevated CO2") +
  xlab("Treatment") + theme_bw(base_size = 18)  + scale_y_continuous(limits = c(0.98, 1.00)) + theme(axis.text.x=element_blank(),
                                                       axis.title.x=element_blank(),
                                                       axis.title.y=element_blank(),
                                                       legend.position="none")
eCO2_CN
ggsave(filename = "eCO2_CN_prim.png",
       plot = eCO2_CN, #this is what you named your plot as
       bg = "transparent",
       width = 2.5, height = 2, units = "in",
       dpi = 300)

#NPP+C:N - site clay and both treatments
level_order <- c('Bulk_CtoN', 'MicC', 'CopiotoOligo', 'SOMc','SOMp', 'N_min')
eCO2_pools.trt <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(SOMp = SOM_1) %>%
  mutate(SOMc = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(CopiotoOligo = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, MicC, CopiotoOligo, SOMc, SOMp, N_min) %>%
  pivot_longer(cols= 4:9, values_to = "dif", names_to = "Pool") %>%
  filter(tx != "Baseline") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, color = tx, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) + ylim(0.8,1.3) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools", color='Treatment', shape="Allocation") +
  theme_bw(base_size=16)
eCO2_pools.trt
#NPP+C:N - site clay and solely priming
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
#NPP - site clay and solely priming
level_order <- c('Bulk_CtoN', 'MicC', 'Micr_K', 'SOMc', 'N_min')
eCO2_pools.prim2 <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(SOMp = SOM_1) %>%
  mutate(SOMc = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(Micr_K = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, MicC, Micr_K, SOMc, N_min) %>%
  pivot_longer(cols= 4:8, values_to = "dif", names_to = "Pool") %>%
  filter(tx == "Priming") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) + ylim(0.8,1.1) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools", shape="Allocation to \nroot exudates") +
  theme_bw(base_size=16)
eCO2_pools.prim2
#NPP+C:N - 5% clay
level_order <- c('Bulk_CtoN', 'MicC', 'CopiotoOligo', 'SOMc','SOMp', 'N_min')
eCO2_pools.5 <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(SOMp = SOM_1) %>%
  mutate(SOMc = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(CopiotoOligo = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, MicC, CopiotoOligo, SOMc, SOMp, N_min) %>%
  pivot_longer(cols= 4:9, values_to = "dif", names_to = "Pool") %>%
  filter(tx != "Baseline") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, color = tx, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) + ylim(0.8,1.3) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools", color='Treatment', shape="Allocation") +
  geom_text(x=1.25, y=1.3, label="a) 5% clay", color="black", size=5) + theme_bw(base_size=16)
eCO2_pools.5
#NPP+C:N - 55% clay
level_order <- c('Bulk_CtoN', 'MicC', 'CopiotoOligo', 'SOMc','SOMp', 'N_min')
eCO2_pools.55 <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(SOMp = SOM_1) %>%
  mutate(SOMc = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(CopiotoOligo = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, MicC, CopiotoOligo, SOMc, SOMp, N_min) %>%
  pivot_longer(cols= 4:9, values_to = "dif", names_to = "Pool") %>%
  filter(tx != "Baseline") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, color = tx, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) + ylim(0.8,1.3) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools", color='Treatment', shape="Allocation") +
  geom_text(x=1.25, y=1.3, label="b) 55% clay", color="black", size=5) + theme_bw(base_size=16)
eCO2_pools.55
ggsave(filename = "eCO2_prim_noNPP.png",
       plot = eCO2_pools.prim2, #this is what you named your plot as
       bg = "transparent",
       width = 8.2, height = 5, units = "in",
       dpi = 300)
ggsave(filename = "eCO2_clay5_50_check.png",
       plot = eCO2_pools.5, #this is what you named your plot as
       bg = "transparent",
       width = 8.2, height = 5, units = "in",
       dpi = 300)
ggsave(filename = "eCO2_clay55_50_check.png",
       plot = eCO2_pools.55, #this is what you named your plot as
       bg = "transparent",
       width = 8.2, height = 5, units = "in",
       dpi = 300)

#NPP
level_order <- c('Bulk_CtoN', 'Oligotrophs', 'CopiotoOligo', 'POM','MAOM', 'N_min')
eCO2_pools.npp <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(MAOM = SOM_1) %>%
  mutate(POM = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(CopiotoOligo = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, Oligotrophs, CopiotoOligo, POM, MAOM, N_min) %>%
  pivot_longer(cols= 4:9, values_to = "dif", names_to = "Pool") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, color = tx, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) +
  scale_color_manual(values=c("#009E73", "#56B4E9","#E69F00")) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools") +
  geom_text(x=1.2, y=0.8, label="b) +20% NPP", color="black") + theme_bw(base_size=16)
eCO2_pools.npp
#litter C:N
level_order <- c('Bulk_CtoN', 'Oligotrophs', 'CopiotoOligo', 'POM', 'MAOM', 'N_min')
eCO2_pools.CN <- df_diff %>%
  mutate(fMAOM = SOM_1/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(fPOM = SOM_2/(SOM_1+SOM_2+SOM_3)) %>%
  mutate(MAOM = SOM_1) %>%
  mutate(POM = SOM_2) %>%
  mutate(Oligotrophs = MIC_2) %>%
  mutate(CopiotoOligo = MicC_O) %>%
  mutate(Bulk_CtoN = bulkCN) %>%
  mutate(N_min = Nmin) %>%
  select(experiment, tx, alloc, Bulk_CtoN, Oligotrophs, CopiotoOligo, POM, MAOM, N_min) %>%
  pivot_longer(cols= 4:9, values_to = "dif", names_to = "Pool") %>%
  ggplot(aes(x = factor(Pool, level=level_order), y = dif, color = tx, group = alloc))  +
  geom_hline(yintercept = 1, color = "black", linewidth=1) +
  geom_point(aes(shape=alloc),size=5, alpha=0.7,position = position_jitter(width=0.25) ) +
  scale_color_manual(values=c("#009E73", "#56B4E9","#E69F00")) +
  labs(y = "Change in pools after \n50 years of elevated CO2", x= " Pools") +
  geom_text(x=1.2, y=0.8, label="c) +10% Litter C:N", color="black") + theme_bw(base_size=16)
eCO2_pools.CN

ggsave(filename = "eCO2_pools.png",
       plot = eCO2_pools, #this is what you named your plot as
       bg = "transparent",
       width = 8, height = 4, units = "in",
       dpi = 300)

