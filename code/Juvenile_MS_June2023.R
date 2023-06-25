# JUV SOSP MANUSCRIPT # 

# Load packages ####

library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)
library(kableExtra)
library(Rmisc)
library(lme4)
library(car)
library(lmtest)
library(lmerTest)

# Load data ####

data <- read.csv("~/Desktop/SG_juv_MS/behaviour/sg_ms_behaviour.csv")
data %>% head()

brain_data <- read_csv("~/Desktop/SG_juv_MS/brain_blood/output_sofia_long_removed_wb.csv")
brain_data %>% head()

raw_data <- read.csv("~/Desktop/SG_juv_MS/peripheral_tissue/03-26-2023_juv_tissue_imputed.csv") %>% select(-tissue_weight, -conc)
raw_data %>% head()

## Data transformation ####

brain_data$log_con <- log(brain_data$concentration)
brain_data %>% summary()
brain_data %>% view()

raw_data$log_con <- log(raw_data$raw) 
raw_data %>% summary()

# **********************************************####

# Behaviour Analyses ####

## Call Latency ####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(call_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(call_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(call_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Song Latency ####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(song_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(song_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(song_lat ~ age, data = ., 
              na.rm = T)

## Flight Latency ####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flight_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flight_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(flight_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Response Latency #####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(resp_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(resp_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(resp_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Number of Songs #####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(songs ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(songs ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(songs ~ age, data = ., 
              na.rm = T, conf.int = T)

## Number of Flights ####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flights ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flights ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(flights ~ age, data = ., 
              na.rm = T, conf.int = T)

## Time Within 5m ####

#### STI juveniles vs controls ####
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(time ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(time ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
data %>%
  filter(condition == "control") %>% 
  wilcox.test(time ~ age, data = ., 
              na.rm = T, conf.int = T)

# **********************************************####

# Steroid Analyses ####

## Testosterone ####

#### Blood ####

mw_t_circ <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_t_circ

#### Brain ####

t.brain.out <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah", "bnst", "cb", "cg", 
                       "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  lmer(log_con ~ matrix*condition + (1|subject), data = .)
  anova(t.brain.out)

t.brain.emm <- 
  emmeans::emmeans(t.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
t.brain.emm

#### Peripheral Tissues ####

tt.out <- raw_data %>% 
  filter(steroid == "t") %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
  anova(tt.out)

t.peripheral.emm <- 
  emmeans::emmeans(tt.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
t.peripheral.emm

## Androstenedione ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

ae.adrenal.out <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
ae.adrenal.out

## Dihydrotestosterone ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

## Dehydroepiandrosterone ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

## 17b-estradiol ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

e2.liver.out <- raw_data %>% 
  filter(steroid == "e2") %>% 
  filter(matrix == "liver") %>% 
  wilcox.test(raw ~ condition, data = .)
e2.liver.out

## 17a-estradiol ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

ae2.adrenal.out <- raw_data %>% 
  filter(steroid == "ae2") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
ae2.adrenal.out

## Estrone ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

## Estriol ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

e3.adrenal.out <- raw_data %>% 
  filter(steroid == "e3") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
e3.adrenal.out

## Progesterone ####

#### Blood ####

mw_p4_circ <- brain_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_p4_circ

#### Brain ####

#### Peripheral Tissues ####

p4.out <- raw_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(p4.out)

p4.peripheral.emm <- 
  emmeans::emmeans(p4.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
p4.peripheral.emm

## Pregnenolone ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

p5.adrenal.out <- raw_data %>% 
  filter(steroid == "p5") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
p5.adrenal.out

## Corticosterone ####

#### Blood ####


#### Brain ####

b.brain.out <- brain_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls",
                       "ncm", "poa", "tna", "vmh", "vta")) %>% 
  lmer(log_con ~ matrix*condition + (1|subject), data = .)
anova(b.brain.out) 

b.brain.emm <- 
  emmeans::emmeans(b.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite") 
b.brain.emm

#### Peripheral Tissues ####

b.out <- raw_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(b.out)

b.peripheral.emm <- 
  emmeans::emmeans(b.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
b.peripheral.emm

## Dehydrocorticosterone ####

#### Blood ####

mw_dhc_circ <- brain_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_dhc_circ  

#### Brain ####

dhc.brain.out <- brain_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  lmer(log_con ~ condition * matrix + (1|subject), data = .)
  anova(dhc.brain.out)

dhc.brain.emm <- 
  emmeans::emmeans(dhc.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
dhc.brain.emm

#### Peripheral Tissues ####

dhc.out <- raw_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(dhc.out)

dhc.peripheral.emm <- 
  emmeans::emmeans(dhc.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")
dhc.peripheral.emm

# **********************************************####






