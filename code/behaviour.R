# BEHAVIOUR 
#copied from 2-21-2022 Juv beh data (barplots) file
#21 Dec 2022

library(tidyverse)
library(rstatix)
library(ggpubr) 
library(kableExtra)
library(Rmisc)
library(lme4) ; library(car)
library(lmtest)

## load data ####  
path <- c("SG_juv_MS/behaviour/")
data <- read.csv(paste0(path, "sg_ms_behaviour.csv")) 
data %>% head()

## shapiro test 
data %>% dplyr::select(call_lat:time) %>% 
  sapply(., shapiro.test)

## data visualization ####
data %>% select(call_lat:time) %>% 
  gather(key, value) %>% 
  mutate(ln_val = log1p(value)) %>% 
  ggplot(aes(x = ln_val)) + geom_histogram() +
  theme_classic() +
  facet_wrap(~key, scales = "free")

data %>% select(condition:time) %>% 
  gather(key, value, -condition, - age) %>% 
  mutate(ln_val = log1p(value)) %>% 
  ggplot(aes(x = condition, y = ln_val)) +
  geom_boxplot() + theme_classic() +
  facet_wrap(~key, scales = "free")

data %>% select(condition:time) %>% 
  gather(key, value, -condition, - age) %>% 
  mutate(ln_val = log1p(value)) %>% 
  ggplot(aes(x = age, y = ln_val)) +
  geom_boxplot() + theme_classic() +
  facet_wrap(~key, scales = "free")

data %>% select(condition:time) %>% 
  gather(key, value, -condition, - age) %>% 
  mutate(ln_val = log1p(value)) %>% 
  ggplot(aes(x = age, y = value, fill = condition)) +
  geom_boxplot() + theme_classic() +
  facet_wrap(~key, scales = "free")

data %>% head()

# Non-parametric analysis #####
#note: adult controls - call latency missing for those that were noted to have called
#note: therefore inappropriate to include control adults in latency comparisons

## call latency ####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(call_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(call_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## song latency ####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(song_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(song_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## flight latency ####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flight_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flight_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## response latency #####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(resp_lat ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(resp_lat ~ age, data = ., 
              na.rm = T, conf.int = T)

## number of songs #####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(songs ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(songs ~ age, data = ., 
              na.rm = T, conf.int = T)

## number of flights ####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(flights ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(flights ~ age, data = ., 
              na.rm = T, conf.int = T)

## time within 5M ####

# mann-whitney test comparing STI juveniles vs controls
data %>% 
  filter(age == "juv") %>% 
  wilcox.test(time ~ condition, data = ., 
              na.rm = T, conf.int = T)

# mann-whitney test comparing sti juveniles vs adults
data %>%
  filter(condition == "sti") %>% 
  wilcox.test(time ~ age, data = ., 
              na.rm = T, conf.int = T)

# ------------------------ #
# Other models: linear models and two-way anovas on log-tranformed data ####
## call latency ####
model1 <- lm((call_lat)~condition*age, data = data)
summary(model1)

aov1 <- aov((call_lat)~condition*age, data = data)
summary(aov1)

#log transformed data
aov1_ln <- aov(log1p(call_lat)~condition*age, data = data)
summary(aov1_ln)
#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(call_lat) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(call_lat) ~ age, data = .)

## song latency ####
model2 <- lm((song_lat)~condition*age, data = data)
summary(model2)

aov2 <- aov((song_lat)~condition*age, data = data)
summary(aov2)

#log transformed data
aov2_ln <- aov(log1p(song_lat)~condition*age, data = data)
summary(aov2_ln)
#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(song_lat) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(song_lat) ~ age, data = .)

## flight latency ####
model3 <- lm((flight_lat)~condition*age, data = data)
summary(model3)

aov3 <- aov((flight_lat)~condition*age, data = data)
summary(aov3)

#log transformed data
aov3_ln<- aov(log1p(flight_lat)~condition*age, data = data)
summary(aov3_ln)

#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(flight_lat) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(flight_lat) ~ age, data = .)

## response latency #####
model4 <- lm((resp_lat)~condition*age, data = data)
summary(model4)

aov4 <- aov((resp_lat)~condition*age, data = data)
summary(aov4)

#log transformed data
aov4_ln <- aov(log1p(resp_lat)~condition*age, data = data)
summary(aov4_ln)

#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(resp_lat) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(resp_lat) ~ age, data = .)

## number of songs #####
model5 <- lm((songs)~condition*age, data = data)
summary(model5)

aov5 <- aov((songs)~condition*age, data = data)
summary(aov5)

#log transformed data
aov5_ln <- aov(log1p(songs)~condition*age, data = data)
summary(aov5_ln)

#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(songs) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(songs) ~ age, data = .)

## number of flights ####
model6 <- lm((flights)~condition*age, data = data)
summary(model6)

aov6 <- aov((flights)~condition*age, data = data)
summary(aov6)

#log transformed data
aov6_ln <- aov(log1p(flights)~condition*age, data = data)
summary(aov6_ln)


#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(flights) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(flights) ~ age, data = .)

## time within 5M ####

model7 <- lm((time)~condition*age, data = data)
summary(model7)

aov7 <- aov((time)~condition*age, data = data)
summary(aov7)
#log transformed data
aov7_ln <- aov(log1p(time)~condition*age, data = data)
summary(aov7_ln)

#t test comparing STI juveniles vs controls
data %>% filter(age=="juv") %>% 
  t.test(log1p(time) ~ condition, data = .)

#t test comparing STI adults vs juveniles
data %>% filter(condition=="sti") %>% 
  t.test(log1p(time) ~ age, data = .)
