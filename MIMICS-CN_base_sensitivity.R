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
# MIMICS single point function - Figure 3b
###########################################
experiment = c('Control',
               'CN=37.5,LIG=15',
               'CN=62.5,LIG=25',
               'CN=25,LIG=10',
               'CN=74.9,LIG=30')
for (e in 1:5)  {   # loop over litter quality
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

  if (e==2) { CN=37.5
     LIG=15
  }
  if (e==3) { CN=62.5
     LIG=25
  }
  if (e==4) { CN=25
    LIG=10
  }
  if (e==5) {CN=74.9
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

CN_LQ.1 <- df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  mutate(CN = c(50, 37.5, 62.5, 25, 74.9)) %>%
  mutate(Lignin = c(20,15,25,10,30)) %>%
  ggplot(aes(x=CN,
             y=CNRatio,show.legend = FALSE)) + #color=experiment,
  geom_point(size = 4) +
  #geom_point(aes(x=Lignin, y=CNratio), size = 4) +
  coord_cartesian(ylim=c(8,10), xlim = c(20,80)) +
  scale_x_continuous(sec.axis=sec_axis(~., breaks=c(50, 37.5, 62.5, 25, 74.9), labels=c('20','15','25','10','30'))) +
  #scale_color_manual(values=c("#009E73", "#56B4E9","#E69F00", "#F0E442","#CC79A7")) +
  labs(y = "Bulk C:N", x = "Litter C:N") +theme_bw(base_size = 20)
CN_LQ.1


###########################################
# MIMICS single point function - Figure 5b
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

