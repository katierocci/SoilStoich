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

#bring in RXEQ function
source("CN_RXEQ.R")
source("calc_Tpars.R")


###########################################
# MIMICS single point function
###########################################
experiment = c('Control',
               'Clay=5','Clay=55',
               'CN=25,LIG=10',
               'CN=75,LIG=30')
for (e in 1:5)  {   # loop over exudation experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ

  Site = 'Temperate Forest'
  ANPP = df$ANPP
  TSOI = df$MAT
  CLAY = df$CLAY2
  LIG = df$LIG
  CN = df$CN
  exud = 0.
  x=1
  ############################################################
  # select range of values to plot
  ############################################################
  CLAY=25
  CN=50
  LIG=20

  if (e==2) {CLAY=5}
  if (e==3) {CLAY=55}
  if (e==4) { CN=25
    LIG=10
  }
  if (e==5) {CN=75
    LIG=30
  }

  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud) #>> Same example site input as used for stode
  Tpars

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
  if (e == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[e]
  } else {
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[e]
    df_ss <- rbind(df_ss, df_exp)
  }

}  # close experiment loop

#====================================
## Finished steady state calculations
## plot results
#====================================

df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (gC/m2)",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

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
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

#replot with just litter changes and aesthetic to match other figures
exp_keep = c('CN=25,LIG=10', 'Control', 'CN=75,LIG=30')
tiff("MIMICS_LQ.tiff", units="px", width=2200, height=2000, res=300)
CN_LQ <- df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  filter(experiment %in% exp_keep) %>%
  ggplot(aes(x=factor(experiment, level = c('CN=25,LIG=10', 'Control', 'CN=75,LIG=30')), y=CNRatio,fill=experiment,show.legend = FALSE)) +
  geom_bar(stat = "identity",show.legend = FALSE) +
  coord_cartesian(ylim=c(6,12)) + scale_fill_manual(values=c("#009E73", "#56B4E9","#E69F00")) +
  labs(y = "Bulk C:N") +theme_bw(base_size = 16) +theme(axis.title.x=element_blank())
CN_LQ
dev.off()
ggsave(filename = "CN_LQ.png",
       plot = CN_LQ, #this is what you named your plot as
       bg = "transparent",
       width = 6, height = 4, units = "in",
       dpi = 300)

###########################################
# MIMICS single point function
###########################################
experiment = c('high clay, high quality',
               'control',
               'low clay, low quality')
for (e in 1:3)  {   # loop over experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ

  Site = 'Temperate Forest'
  ANPP = df$ANPP
  TSOI = df$MAT
  CLAY = df$CLAY2
  LIG = df$LIG
  CN = df$CN
  exud = 0.
  x=1
  ############################################################
  # select range of values to plot
  ############################################################
  CLAY=25
  CN=50
  LIG=20

  if (e==1) {
    CLAY=55
    CN=25
    LIG=10
  }
  if (e==3) {
    CLAY=5
    CN=75
    LIG=30
  }

  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud) #>> Same example site input as used for stode
  Tpars

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
  if (e == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[e]
  } else {
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[e]
    df_ss <- rbind(df_ss, df_exp)
  }

}  # close experiment loop

#====================================
## Finished steady state calculations
## plot results
#====================================

df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (gC/m2)",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

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
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

###########################################
# MIMICS single point function - recreating Kats soil texture and MAOM fig
###########################################
experiment = c('Clay=10','Clay=35', 'Clay=60', 'Clay=85')
for (e in 1:4)  {   # loop over exudation experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ

  Site = 'Temperate Forest'
  ANPP = df$ANPP
  TSOI = df$MAT
  CLAY = df$CLAY2
  LIG = df$LIG
  CN = df$CN
  exud = 0.
  x=1
  ############################################################
  # select range of values to plot
  ############################################################
  CLAY=10

  if (e==2) {CLAY=35}
  if (e==3) {CLAY=60}
  if (e==4) {CLAY=85}


  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud) #>> Same example site input as used for stode
  Tpars

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
  if (e == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[e]
  } else {
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[e]
    df_ss <- rbind(df_ss, df_exp)
  }

}  # close experiment loop

#====================================
## Finished steady state calculations
## plot results
#====================================

df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (gC/m2)",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

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
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

#replot to match Kat's figure
df_wide <- pivot_wider(df_ss, values_from = "y", names_from = pools)
CN_tex <- df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  inner_join(df_wide, by = "experiment") %>%
  mutate(fMAOM = (SOM_p/(SOM_p+SOM_c+SOM_a))*100) %>%
  ggplot(aes(x=experiment, y=CNRatio, show.legend = FALSE)) +
  geom_point(color="black", size =3) + geom_point(aes(x=experiment, y=fMAOM/6.5), color="green3", size =3) +
  scale_y_continuous(name="Bulk soil C:N", limits = c(4,16), sec.axis = sec_axis(~.*5.5, name="MAOM-C/SOM-C (%)"))+
  theme_bw(base_size = 12) +theme(axis.title.x=element_blank(),
                                  axis.title.y.right = element_text(color = "green3"), axis.line.y.right =element_line(color = "green3"),
                                  axis.ticks.y.right =element_line(color = "green3"), axis.text.y.right = element_text(color = "green3"),
                                  axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black"))
CN_tex
ggsave(filename = "CN_tex.png",
       plot = CN_tex, #this is what you named your plot as
       bg = "transparent",
       width = 6, height = 4, units = "in",
       dpi = 300)

