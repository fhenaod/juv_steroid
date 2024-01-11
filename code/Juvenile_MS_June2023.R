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

# behaviour data
beh_data <- read.csv("data/sg_ms_behaviour.csv")
beh_data %>% head()

# brain and circulation data
brain_data <- read_csv("data/sg_ms_brain.csv")
brain_data %>% head()

# peripheral tissue data 
pt_data <- read.csv("data/sg_ms_periph_tissue.csv") %>% 
  select(-tissue_weight)
pt_data %>% head()

## Data transformation ####
brain_data$log_con <- log(brain_data$concentration)
brain_data %>% summary()
brain_data %>% view()

pt_data$log_con <- log(pt_data$conc) 
pt_data %>% summary()

# **********************************************####

# Behaviour Analyses ####

## Call Latency ####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(call_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(call_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(call_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Song Latency ####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(song_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(song_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(song_lat ~ age, data = ., 
              na.rm = T)

## Flight Latency ####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flight_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flight_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(flight_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Response Latency #####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(resp_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(resp_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(resp_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## Number of Songs #####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(songs ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(songs ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(songs ~ age, data = ., 
              na.rm = T, conf.int = T)

## Number of Flights ####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flights ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flights ~ age, data = ., 
              na.rm = T, conf.int = T)

#### Control juveniles vs adults ####
beh_data %>%
  filter(condition == "control") %>% 
  wilcox.test(flights ~ age, data = ., 
              na.rm = T, conf.int = T)

## Time Within 5m ####

#### STI juveniles vs controls ####
beh_data %>% 
  filter(age == "juv") %>% 
  wilcox.test(time ~ condition, data = ., 
              na.rm = T, conf.int = T)

#### STI juveniles vs adults ####
beh_data %>%
  filter(condition == "sti") %>% 
  wilcox.test(time ~ age, data = ., 
              na.rm = T, conf.int = T)

#### rater's differences ####
library(broom)
beh_data %>% 
  select_if(is.numeric) %>% 
  sapply(., shapiro.test)

beh_data %>% dplyr::select(rater = Rater, call_lat:time) %>% 
  gather(key, value, -rater) %>% 
  group_by(key) %>% 
  ggplot(aes(x = rater, y = value)) + 
  geom_boxplot(notchwidth = T) + 
  geom_point(col = "red") +
  theme_minimal() + stat_compare_means() +
  facet_wrap(~key, scales ="free")

#### Control juveniles vs adults ####
beh_data %>%
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

tt.out <- pt_data %>% 
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

ae.adrenal.out <- pt_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(conc ~ condition, data = .)
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

e2.liver.out <- pt_data %>% 
  filter(steroid == "e2") %>% 
  filter(matrix == "liver") %>% 
  wilcox.test(conc ~ condition, data = .)
e2.liver.out

## 17a-estradiol ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

ae2.adrenal.out <- pt_data %>% 
  filter(steroid == "ae2") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(conc ~ condition, data = .)
ae2.adrenal.out

## Estrone ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

## Estriol ####

#### Blood ####

#### Brain ####

#### Peripheral Tissues ####

e3.adrenal.out <- pt_data %>% 
  filter(steroid == "e3") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(conc ~ condition, data = .)
e3.adrenal.out

## Progesterone ####

#### Blood ####

mw_p4_circ <- brain_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_p4_circ

mw_p4_plasma <- brain_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("plasma")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_p4_plasma


#### Brain ####

#### Peripheral Tissues ####

p4.out <- pt_data %>% 
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

p5.adrenal.out <- pt_data %>% 
  filter(steroid == "p5") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(conc ~ condition, data = .)
p5.adrenal.out

## Corticosterone ####

#### Blood ####
mw_b_circ <- brain_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix == "wb") %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_b_circ 

mw_b_plasma <- brain_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix == "wb") %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_b_plasma

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

b.out <- pt_data %>% 
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

mw_dhc_plasma <- brain_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("plasma")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_dhc_plasma

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

dhc.out <- pt_data %>% 
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

