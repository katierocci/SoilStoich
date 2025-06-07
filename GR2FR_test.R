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

#########################################
#steady state calculations for grassland pools
###################################

data = read.csv("LTER_SITE_1.csv")
df <- data[5,] #6=Cedar Creek

Site = df$Site
ANPP = df$ANPP
TSOI = df$MAT
CLAY = df$CLAY2 #5 or 55
LIG = df$LIG
CN = df$CN #45
x=1 #see calc_Tpars: x=1 normal, x=2 exudation x=3 exudation and desorption, x=4 ECM N mining
exud = 0.1

############################################################
# MIMICS MODEL CODE STARTS HERE
############################################################


Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                   x=x,exud=exud,nUPmod=1)


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

df_ss = as.data.frame(test) * MICROtoECO
df_ss = t(df_ss)
colnames(df_ss) = pools
df_gr = as.data.frame(df_ss)

###########################################
# steady state calculations of forest from grassland values
###########################################

  data = read.csv("LTER_SITE_1.csv")
  df <- data[5,] #Cedar Creek but with Andrews Forest litter quantity and quality
  df[,5] <- data[7,5] #ANPP
  df[,11] <- data[7,11] #lignin
  df[,12] <- data[7,12] #N
  df[,13] <- data[7,13] #C:N
  df[,14] <- data[7,14] #relANPP


  Site = df$Site
  ANPP = df$ANPP
  TSOI = df$MAT
  CLAY = df$CLAY2 #5 or 55
  LIG = df$LIG
  CN = df$CN #45
  x=1 #see calc_Tpars: x=1 normal, x=2 exudation x=3 exudation and desorption, x=4 ECM N mining
  exud = 0.3

  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################


  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud,nUPmod=1)


  #----------initialize pools---------------
  LIT_1  <<- df_gr$LIT_m
  LIT_2  <<- df_gr$LIT_s
  MIC_1  <<- df_gr$MIC_r
  MIC_2  <<- df_gr$MIC_K
  SOM_1  <<- df_gr$SOM_p
  SOM_2  <<- df_gr$SOM_c
  SOM_3  <<- df_gr$SOM_a

  LIT_1_N  <<- df_gr$LIT_m_N
  LIT_2_N  <<- df_gr$LIT_s_N
  MIC_1_N  <<- df_gr$MIC_r_N
  MIC_2_N  <<- df_gr$MIC_K_N
  SOM_1_N  <<- df_gr$SOM_p_N
  SOM_2_N  <<- df_gr$SOM_c_N
  SOM_3_N  <<- df_gr$SOM_a_N
  DIN      <<- df_gr$DIN
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

  df_ss4 = as.data.frame(test) * MICROtoECO
  df_ss4 = t(df_ss4)
  colnames(df_ss4) = pools
  df_fr_nomining = as.data.frame(df_ss4)

  #bringing together for plotting
  df_gr$trt <- "grassland"
  df_fr_old$trt <- "coniferous forest"
  df_fr$trt <- "grassland_to_forest"
  df_fr_nomining$trt <- "G2F_nomining"
  df_all <- rbind(df_gr, df_fr_old, df_fr, df_fr_nomining)


  df_all %>%
    pivot_longer(1:15, names_to = "pools", values_to = "value") %>%
    ggplot(aes(x = pools, y = value, color = trt))  +
    geom_point(size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
    scale_color_manual(values = c("lightskyblue", "darkgoldenrod2", "aquamarine4", "#CC79A7")) +
    labs(y = "carbon content")+theme_bw()

  df_all %>%
    select(SOM_p, SOM_c, SOM_a, SOM_p_N, SOM_c_N, SOM_a_N, trt) %>%
    pivot_longer(1:6, names_to = "pools", values_to = "value") %>%
    ggplot(aes(x = pools, y = value, color = trt))  +
    geom_point(size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
    scale_color_manual(values = c("lightskyblue", "darkgoldenrod2", "aquamarine4", "#CC79A7")) +
    labs(y = "Carbon conent")+theme_bw()

  df_all %>%
    select(SOM_p, SOM_c, SOM_a, SOM_p_N, SOM_c_N, SOM_a_N, trt) %>%
    filter(trt != "coniferous forest") %>%
    pivot_longer(1:6, names_to = "pools", values_to = "value") %>%
    pivot_wider(values_from = "value", names_from = "trt") %>%
    mutate(RR = grassland_to_forest/grassland) %>%
    ggplot(aes(x = pools, y = RR))  +
    geom_point(size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
    geom_hline(yintercept = 1, color = "black", linewidth = 1) +
    labs(y = "Ratio of C in grassland to forest \n relative to grassland")+theme_bw(base_size = 16)


  tiff("GR2FR.tiff", units="px", width=2800, height=1500, res=300)
  df_all %>%
    select(SOM_p, SOM_c, SOM_a, SOM_p_N, SOM_c_N, SOM_a_N, trt) %>%
    filter(trt != "coniferous forest") %>%
    pivot_longer(1:6, names_to = "pools", values_to = "value") %>%
    pivot_wider(values_from = "value", names_from = "trt") %>%
    mutate(RR1 = grassland_to_forest/grassland) %>%
    mutate(RR2 = G2F_nomining/grassland) %>%
    ggplot()  +
    geom_point(aes(x = pools, y = RR1, color = "With ECM N mining"), size=5, alpha=0.7, position = position_jitter(width=0.1) ) +
    geom_point(aes(x = pools, y = RR2, color = "Without ECM N mining"), size=5, alpha=0.7,position = position_jitter(width=0.1) ) +
    geom_hline(yintercept = 1, color = "black", linewidth = 1) +
    labs(y = "Ratio of C in grassland to forest \n relative to grassland")+theme_bw(base_size = 16) +
    scale_color_manual(name='Model type',
                       breaks=c('With ECM N mining', 'Without ECM N mining'),
                       values=c('With ECM N mining'='aquamarine4', 'Without ECM N mining'='#CC79A7'))
  dev.off()
