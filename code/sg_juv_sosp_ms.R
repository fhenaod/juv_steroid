#juv SOSP manuscript 2022

#GOAL AND IMPORTANT NOTES####
#goal: juvenile song sparrow manuscript figures + analyses for submission March 2023

#desired analyses: 
  # major components: 
  # (1) mixed anova OR lmm analyzing the effect of condition (levels: STI or control) 
    #and matrix (levels: a:10 brain regions, b: wb, plasma, c: liver, pectoral muscle, adrenals, testes) 
    #on steroid concentration, including interaction
    #post-hoc pairwise analyses on significant main effects and interactions with benjamini hochberg (?) correction
    #steroids analyzed separately
  # (2) PCA on behavioural data --> create composite 'aggression score' based on individual coordinates
  # (3) in STI subjects only, create correlation matrix looking at 
    #relationship between 'aggression score' coordinates and steroid concentrations across brain regions
    #steroids analyzed separately
  # supplemental figures/analyses: 
  # (1) correlation matrices of different steroid concentrations across brain regions (i.e. B and DHC)
  # (2) tables of steroid concentrations in peripheral tissues
    #analyzed by (nonparametic - mann whitney if necessary) t-test if steroid is detectable in one tissue only
    #analysed by mixed anova OR lmm if steroid is detectable in multiple tissues
    #analyze steroids separately; group tables by steroid class (androgens, estrognens, progestins, glucocorticoids)

#packages
library(tidyverse)
library(rstatix)
library(ggpubr) ; library(patchwork)
library(kableExtra)
library(Rmisc)
library(lme4) ; library(car)
library(lmtest)


#BRAIN AND BLOOD ANALYSES####

#NOV 2022 based on ng/g or ng/ml concentrations with imputation threshold set to 20% detectable
#CJ124 (RT shift) and JUV 1 TNA + VTA (extraction error) removed from output_sofia_long --> now called output_sofia_long_removed
  #for CJ124, these matrices have to be removed (this is done for AE, T, B, and DHC in sofia_output_long_removed.csv): 
    #AE:    tna, ncm, vta, cg, cb, wb, plasma
    #DHEA:  ncm, vta, cg, cb, wb, plasma
    #T:     ncm, vta, cg, cb, wb plasma
    #B:     vta, cg, cb, wb, plasma
    #DHC:   cb, wb, plasma
    #E2:    vta, cg, cb, wb, plasma
    #aE2:   vta, cg, cb, wb, plasma
    #E1:    vta, cg, cb, wb, plasma

#removing NAC from everything, 21 nov 2022 - we will not be showing NAC data in MS

#1 Dec 2022: In a new csv file (output_sofia_long_removed_wb), I manually converted WB concentration from ng/ml (ng/.005ml) which was used for imputation 
  #to ng/g (ng/.0051145g) based off of Kiran's recommendation and new WB weight based on Jerome measurement (1ml = 1.0229g)
    #I multiplied wb concentration column by .005 (gets back to just ng), then divided ng by .0051145 (this is how much 5ul weighs according to Jerome)

# load data ####
path <- c("SG_juv_MS/brain_blood/")
brain_data <- read_csv(paste0(path, "output_sofia_long_removed_wb.csv"))

path4figs <- c("SG_juv_MS/r_figs/")

# Data visualization, exploration and transformation ####
## qqplots and histograms 
t_brain <- brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
qqnorm(t_brain$concentration) ; qqline(t_brain$concentration, col = "red")
hist(t_brain$concentration)

t_blood <- brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("wb", "plasma"))
qqnorm(t_blood$concentration) ; qqline(t_blood$concentration, col = "red")
hist(t_blood$concentration)


b_brain <- brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
qqnorm(b_brain$concentration); qqline(b_brain$concentration, col = "red")
hist(b_brain$concentration)

b_blood <- brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("wb", "plasma"))
qqnorm(b_blood$concentration); qqline(b_blood$concentration, col = "red")
hist(b_blood$concentration)


dhc_brain <- brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
qqnorm(dhc_brain$concentration) ; qqline(dhc_brain$concentration, col = "red")
hist(dhc_brain$concentration)

dhc_blood <- brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb", "plasma"))
qqnorm(dhc_blood$concentration) ; qqline(dhc_blood$concentration, col = "red")
hist(dhc_blood$concentration)


p4_brain <- brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
qqnorm(p4_brain$concentration) ; qqline(p4_brain$concentration, col = "red")
hist(p4_brain$concentration)

p4_blood <- brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb", "plasma"))
qqnorm(p4_blood$concentration) ; qqline(p4_blood$concentration, col = "red")
hist(p4_blood$concentration)

## shapiro wilk test of normality (by neural and circulating steroid groups to be analyzed separately)
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(concentration) %>% sapply(., shapiro.test)
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(concentration) %>% sapply(., shapiro.test)

brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(concentration) %>% sapply(., shapiro.test)
brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(concentration) %>% sapply(., shapiro.test)

brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(concentration) %>% sapply(., shapiro.test)
brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(concentration) %>% sapply(., shapiro.test)

brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(concentration) %>% sapply(., shapiro.test)
brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(concentration) %>% sapply(., shapiro.test)

## levene's test for homogeneity of variance
    #do you include matrix?? warning message about group coercion without it
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .) #good
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .) #good

brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .)
brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .)

brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .) #good
brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .) #good

brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .)
brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(concentration ~ condition*matrix, data = .)

## add new column of log transformed data 
# since data are not normal nor homogeneous (see * exceptions above) 
brain_data$log_con <- log(brain_data$concentration)
brain_data %>% summary()

#re-run shapiro wilk test of normality after log transformation
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #still not normal
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #still not normal

brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #much improved, but still not normal: p = .01
brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #normal after transformation

brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #still not normal
brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #still not normal

brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #still not normal
brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  select(log_con) %>% sapply(., shapiro.test) #much improved, but still not normal: p = .003

#re-run levene's test for homogeneity of variance 
#do you include matrix?? warning message about group coercion without it
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #got worse - now not homogeneous
brain_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #still good

brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #good now
brain_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #good now

brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #still good
brain_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #still good

brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #still not homogeneous
brain_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("wb", "plasma")) %>% 
  leveneTest(log_con ~ condition*matrix, data = .) #good now

# Analyses ####

#do NOT run stats with NAs as zero

## testosterone ####
# anova is throwing error b/c AH, VMH, VTA, BNST, and Cb is detectable ONLY in STI or control group - should I exclude from analyses then?  
res_aov_t <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_t)

#see - if I remove these, then the warning goes away...
res_aov_t <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("cg", "ls", "ncm", "poa", "tna")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_t)

#I know I shouldn't do these b/c no significant effects from ANOVA table- just to set up code
#pairwise posthoc: sti vs control grouped by matrix
pwc_t <-  brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_t

#lmm: also can't handle that some matrices have NA for an entire condition; lmm auto-drops these columns
#here, interaction is significant (but sketchy); F value not that different btw ANOVA and LMM - maybe correction is changing p value
t <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah", "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
t.out <- (lmer(log_con ~ matrix*condition + (1|subject), data = t))
anova(t.out)
hist(residuals(t.out))
shapiro.test(residuals(t.out))
qqnorm(residuals(t.out))

#mann whitney on blood testosterone
mw_t_circ <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_t_circ

## corticosterone ####
res_aov_b <- brain_data %>% 
    filter(steroid == "b") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
    anova_test(data = ., dv = log_con, wid = subject,
               within = matrix, between = condition)
get_anova_table(res_aov_b)
 
#pairwise posthoc: sti vs control grouped by matrix 
pwc_b <-  brain_data %>% 
    filter(steroid == "b") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%   
    group_by(matrix) %>%
    pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_b
  
#lmm - no warning messages (like with t); matches anova results
b <- brain_data %>% 
    filter(steroid == "b") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
anova(lmer(log_con ~ matrix*condition + (1|subject), data = b)) 
  
#mann whitney on blood corticosterone
mw_b_circ <- brain_data %>% 
    filter(steroid == "b") %>% 
    filter(matrix %in% c("wb")) %>% 
    wilcox.test(concentration ~ condition, data = .)
mw_b_circ  

## dehydrocorticosterone ####
# anova is throwing error b/c Cb is detectable ONLY in STI group - should I exclude from analyses then?  
res_aov_dhc <- brain_data %>% 
    filter(steroid == "dhc") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
    anova_test(data = ., dv = log_con, wid = subject,
               within = matrix, between = condition)
  get_anova_table(res_aov_dhc)
  
  #cb removed; if you decide to remove for the ANOVA, also remove for posthoc (see below)
  res_aov_dhc <- brain_data %>% 
    filter(steroid == "dhc") %>% 
    filter(matrix %in% c("ah" , "bnst", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
    anova_test(data = ., dv = log_con, wid = subject,
               within = matrix, between = condition)
  get_anova_table(res_aov_dhc)
  
  #pairwise posthoc: sti vs control grouped by matrix 
  pwc_dhc <-  brain_data %>% 
    filter(steroid == "dhc") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%   
    group_by(matrix) %>%
    pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
  pwc_dhc
  
  #lmm - warning message about Cb column getting dropped (like with t); 
  #ANOVA found both simple + interaction, but lmm found no effect of condition
  dhc <- brain_data %>% 
    filter(steroid == "dhc") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
  anova(lmer(log_con ~ matrix*condition + (1|subject), data = dhc))

  #mann whitney on blood corticosterone
  mw_dhc_circ <- brain_data %>% 
    filter(steroid == "dhc") %>% 
    filter(matrix %in% c("wb")) %>% 
    wilcox.test(concentration ~ condition, data = .)
  mw_dhc_circ  
  
## progesterone ####
  #only running repeated measures ANOVA on matrix b/c p4 only detectable in STI group (but, do I even care about differences across matrices?)
  res_aov_p4 <- brain_data %>% 
    filter(steroid == "p4") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
    anova_test(data = ., dv = log_con, wid = subject,
               within = matrix)
  get_anova_table(res_aov_p4)
  
  #pairwise posthoc: differences across matrices...again I don't really care - also, it's not working...
  pwc_p4 <-  brain_data %>% 
    filter(steroid == "p4") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
    group_by(condition) %>% 
    pairwise_t_test(log_con ~ matrix, p.adjust.method = "BH")
  pwc_p4
  
  #lmm - again, I dropped condition as a fixed effect b/c p4 only detectable in STI group
  p4 <- brain_data %>% 
    filter(steroid == "p4") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta"))
  anova(lmer(log_con ~ matrix + (1|subject), data = dhc))

#mann whitney on blood progesterone
  mw_p4_circ <- brain_data %>% 
    filter(steroid == "p4") %>% 
    filter(matrix %in% c("wb")) %>% 
    wilcox.test(concentration ~ condition, data = .)
  mw_p4_circ


#********************************************************************************####  
# BEHAVIOUR ####
#copied from 2-21-2022 Juv beh data (barplots) file
#21 Dec 2022

# load data ####  
path <- c("SG_juv_MS/behaviour/")
data <- read.csv(paste0(path, "sg_ms_behaviour.csv")) 
data %>% head()

## shapiro test ####
data %>% select(call_lat:time) %>% 
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

## LINEAR MODELS AND TWO-WAY INTERACTIVE ANOVAS FOR EACH BEHAVIOUR #####
#note: adult controls - call latency missing for those that were noted to have called
#note: therefore inappropriate to include control adults in latency comparisons

#call latency
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

#song latency
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

#flight latency
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

#response latency
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

#number of songs
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

#number of flights
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

#time within 5M
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

# PERIPHERAL TISSUES 
#7 November 2022

#the raw column is ng/g (the other columns are incorrect) 
#imputations were made on ng/g values, NOT pg values

path <- c("SG_juv_MS/peripheral_tissue/")
raw_data <- read.csv(paste0(path,"output_data.csv")) %>% 
  select(-tissue_weight, -conc)

#add new column of log transformed data since data are not normally distributed
raw_data$log_con <- log1p(raw_data$raw) 
raw_data %>% summary()

# TABLES ####
#changing to wide format to make tables
#raw concentration
data_wide_raw <- spread(raw_data, matrix, raw) %>% select(-X)
data_wide_raw

#log transformed concentration
data_wide_raw$log_adrenal <- log1p(data_wide_raw$adrenal)
data_wide_raw$log_liver <- log1p(data_wide_raw$liver)
data_wide_raw$log_pec <- log1p(data_wide_raw$pec)
data_wide_raw$log_testes <- log1p(data_wide_raw$testes)
data_wide_raw %>% summary()

#what the heck was I doing here? 4dec2022 sofia doesn't understand previous sofia...
#oh I think I was trying to set up to make summary stat tables by matrix...mhm
data_wide_raw <- spread(raw_data, steroid, log_con) %>% select(-raw)
adrenal_data <- data_wide_steroid_log %>% filter(matrix == "adrenal")
liver_data <- data_wide_steroid_log %>% filter(matrix == "liver")

adrenal_summary <- data.frame(
  mean = c(mean(adrenal_data$ae, na.rm = T),
           mean(adrenal_data$ae2, na.rm = T), 
           mean(adrenal_data$b, na.rm = T), 
           mean(adrenal_data$dhc, na.rm = T), 
           mean(adrenal_data$dht, na.rm = T), 
           mean(adrenal_data$e2, na.rm = T), 
           mean(adrenal_data$e3, na.rm = T), 
           mean(adrenal_data$p4, na.rm = T),
           mean(adrenal_data$p5, na.rm = T),
           mean(adrenal_data$t, na.rm = T)),
  n = c(length(adrenal_data$ae),
        length(adrenal_data$ae2), 
        length(adrenal_data$b), 
        length(adrenal_data$dhc), 
        length(adrenal_data$dht), 
        length(adrenal_data$e2), 
        length(adrenal_data$e3), 
        length(adrenal_data$p4),
        length(adrenal_data$p5),
        length(adrenal_data$t)),
  sd = c(sd(adrenal_data$ae, na.rm = T),
         sd(adrenal_data$ae2, na.rm = T), 
         sd(adrenal_data$b, na.rm = T), 
         sd(adrenal_data$dhc, na.rm = T), 
         sd(adrenal_data$dht, na.rm = T), 
         sd(adrenal_data$e2, na.rm = T), 
         sd(adrenal_data$e3, na.rm = T), 
         sd(adrenal_data$p4, na.rm = T),
         sd(adrenal_data$p5, na.rm = T),
         sd(adrenal_data$t, na.rm = T)),
  row.names = c('AE','AE2', 'B', 'DHC', 'DHT', 'E2', 'E3', 'P4', 'P5', 'T'))

liver_summary <- data.frame(
  mean = c(mean(liver_data$ae, na.rm = T),
           mean(liver_data$ae2, na.rm = T), 
           mean(liver_data$b, na.rm = T), 
           mean(liver_data$dhc, na.rm = T), 
           mean(liver_data$dht, na.rm = T), 
           mean(liver_data$e2, na.rm = T), 
           mean(liver_data$e3, na.rm = T), 
           mean(liver_data$p4, na.rm = T),
           mean(liver_data$p5, na.rm = T),
           mean(liver_data$t, na.rm = T)),
  n = c(length(liver_data$ae),
        length(liver_data$ae2), 
        length(liver_data$b), 
        length(liver_data$dhc), 
        length(liver_data$dht), 
        length(liver_data$e2), 
        length(liver_data$e3), 
        length(liver_data$p4),
        length(liver_data$p5),
        length(liver_data$t)),
  sd = c(sd(liver_data$ae, na.rm = T),
         sd(liver_data$ae2, na.rm = T), 
         sd(liver_data$b, na.rm = T), 
         sd(liver_data$dhc, na.rm = T), 
         sd(liver_data$dht, na.rm = T), 
         sd(liver_data$e2, na.rm = T), 
         sd(liver_data$e3, na.rm = T), 
         sd(liver_data$p4, na.rm = T),
         sd(liver_data$p5, na.rm = T),
         sd(liver_data$t, na.rm = T)),
  row.names = c('AE','AE2', 'B', 'DHC', 'DHT', 'E2', 'E3', 'P4', 'P5', 'T'))
#end of confusion - 4dec2022

#shapiro wilk test of normality  
raw_data %>% filter(steroid == "ae") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>% 
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "ae2") %>% 
  filter(matrix %in% c("adrenal")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "e2") %>% 
  filter(matrix %in% c("liver")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "e3") %>% 
  filter(matrix %in% c("adrenal")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "dht") %>% 
  filter(matrix %in% c("testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "p5") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)

#levene's test for homogeneity of variance 
#do you include matrix?? warning message about group coercion without it
raw_data %>% filter(steroid == "ae") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>% 
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "ae2") %>% 
  filter(matrix %in% c("adrenal")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "e2") %>% 
  filter(matrix %in% c("liver")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "e3") %>% 
  filter(matrix %in% c("adrenal")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "dht") %>% 
  filter(matrix %in% c("testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .) #dht only detectable in sti testes - no need for stats
raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "p5") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)


#summary tables
summarySE(data = raw_data, 
          measurevar = "raw", 
          groupvars = c("condition", "matrix", "steroid"))

raw_data %>% group_by(steroid, matrix, condition) %>% 
  dplyr::summarize(mean = mean(raw, na.rm = T),
                   sd = sd(raw, na.rm = T)) %>% 
  dplyr::arrange(matrix, steroid)

# androstenedione ####
#ae only detectable in both sti and control adrenal - so just a mann whitney on it
# add mann whitney stats to table
mw_ae <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_ae

ae <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% as.data.frame()

#boxplot
raw_data %>% filter(steroid == "ae") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "AE") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "ae") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "Androstenedione") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/ae_peripheral.pdf", width=6, height=4, dpi=100)

# alpha estradiol ####

#ae2 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_ae2 <- raw_data %>% 
  filter(steroid == "ae2") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = ., exact = FALSE)
mw_ae2

#boxplot
raw_data %>% filter(steroid == "ae2") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "aE2") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "ae2") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "alpha estradiol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .4), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/ae2_peripheral.pdf", width=6, height=4, dpi=100)


# estradiol ####
#e2 mann whitney liver - warning message: can't compute exact p-value with ties 
mw_e2 <- raw_data %>% 
  filter(steroid == "e2") %>% 
  filter(matrix == "liver") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e2

#boxplot
raw_data %>% filter(steroid == "e2") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "E2") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "e2") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "Estradiol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/e2_peripheral.pdf", width=6, height=4, dpi=100)


# estriol ####
#e3 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_e3 <- raw_data %>% 
  filter(steroid == "e3") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e3

#boxplot
raw_data %>% filter(steroid == "e3") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Estriol") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))


#barplot
raw_data %>% filter(steroid == "e3") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "estriol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .2), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/e3_peripheral.pdf", width=6, height=4, dpi=100)


# corticosterone ####
#adrenals removed b/c cort concentrations exceed curve
res_aov_b <- raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_b) #effect of matrix, but not conditions


#pairwise posthoc
pwc_b_peripheral <-  raw_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_b_peripheral

#post-hoc with benjamini hochberg - main effect of sti on liver and pec, but not testes
one.way_b <- raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_b

#boxplot
raw_data %>% filter(steroid == "b") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "corticosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/b_peripheral.pdf", width=6, height=4, dpi=100)


# dehydrocorticosterone ####
#adrenals removed b/c dhc conc exceed curve
res_aov_dhc <- raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_dhc) #effect of matrix, but not condition

#pairwise posthoc
pwc_dhc_peripheral <-  raw_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_dhc_peripheral

#post-hoc with benjamini hochberg - main effect of sti on liver 
one.way_dhc <- raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_dhc

#boxplot
raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "dehydrocorticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))


#barplot
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "dehydrocorticosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/dhc_peripheral.pdf", width=6, height=4, dpi=100)


# dihydrotestosterone ####
#dht only detectable in STI testes - no need for stats

#boxplot
raw_data %>% filter(steroid == "dht") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "dht") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "dihydrotestosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/dht_peripheral.pdf", width=6, height=4, dpi=100)


# testosterone ####
res_aov_t <- raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec", "adrenal", "testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_t)

#pairwise posthoc
pwc_t_peripheral <-  raw_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec","adrenal", "testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_t_peripheral

#post hoc looking at effect of sti on t in each matrix (not comparing across matrices)
one.way_t <- raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec", "adrenal", "testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_t

#boxplot
raw_data %>% filter(steroid == "t") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "t") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "testosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/t_peripheral.pdf", width=6, height=4, dpi=100)

# progesterone ####
#adrenals are removed b/c p4 concentration exceeds curve
res_aov_p4 <- raw_data %>% filter(steroid == "p4") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_p4)

#pairwise posthoc with benjamini hochberg correction - looks at differences across matrices
pwc_p4 <- raw_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_p4

#boxplot
raw_data %>% filter(steroid == "p4") %>%  filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "progesterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "progesterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/p4_peripheral.pdf", width=6, height=4, dpi=100)


# pregnenolone ####

#run mann whitney b/c adrenals detectable in sti and control groups; testes only detectable in sti group
mw_p5 <- raw_data %>% 
  filter(steroid == "p5") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_p5


#boxplot
raw_data %>% filter(steroid == "p5") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "pregnenolone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "p5") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "pregnenolone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 120), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/p5_peripheral.pdf", width=6, height=4, dpi=100)


# PERIPHERAL TISSUES #### 
#7 November 2022

#the raw column is ng/g (the other columns are incorrect) 
#imputations were made on ng/g values, NOT pg values

path <- c("SG_juv_MS/peripheral_tissue/")
raw_data <- read.csv(paste0(path,"output_data.csv")) %>% 
  select(-tissue_weight, -conc)

#add new column of log transformed data since data are not normally distributed
raw_data$log_con <- log1p(raw_data$raw) 
raw_data %>% summary()

# TABLES ####
#changing to wide format to make tables
#raw concentration
data_wide_raw <- spread(raw_data, matrix, raw) %>% select(-X)
data_wide_raw

#log transformed concentration
data_wide_raw$log_adrenal <- log1p(data_wide_raw$adrenal)
data_wide_raw$log_liver <- log1p(data_wide_raw$liver)
data_wide_raw$log_pec <- log1p(data_wide_raw$pec)
data_wide_raw$log_testes <- log1p(data_wide_raw$testes)
data_wide_raw %>% summary()

#what the heck was I doing here? 4dec2022 sofia doesn't understand previous sofia...
#oh I think I was trying to set up to make summary stat tables by matrix...mhm
data_wide_raw <- spread(raw_data, steroid, log_con) %>% select(-raw)
adrenal_data <- data_wide_steroid_log %>% filter(matrix == "adrenal")
liver_data <- data_wide_steroid_log %>% filter(matrix == "liver")

adrenal_summary <- data.frame(
  mean = c(mean(adrenal_data$ae, na.rm = T),
           mean(adrenal_data$ae2, na.rm = T), 
           mean(adrenal_data$b, na.rm = T), 
           mean(adrenal_data$dhc, na.rm = T), 
           mean(adrenal_data$dht, na.rm = T), 
           mean(adrenal_data$e2, na.rm = T), 
           mean(adrenal_data$e3, na.rm = T), 
           mean(adrenal_data$p4, na.rm = T),
           mean(adrenal_data$p5, na.rm = T),
           mean(adrenal_data$t, na.rm = T)),
  n = c(length(adrenal_data$ae),
        length(adrenal_data$ae2), 
        length(adrenal_data$b), 
        length(adrenal_data$dhc), 
        length(adrenal_data$dht), 
        length(adrenal_data$e2), 
        length(adrenal_data$e3), 
        length(adrenal_data$p4),
        length(adrenal_data$p5),
        length(adrenal_data$t)),
  sd = c(sd(adrenal_data$ae, na.rm = T),
         sd(adrenal_data$ae2, na.rm = T), 
         sd(adrenal_data$b, na.rm = T), 
         sd(adrenal_data$dhc, na.rm = T), 
         sd(adrenal_data$dht, na.rm = T), 
         sd(adrenal_data$e2, na.rm = T), 
         sd(adrenal_data$e3, na.rm = T), 
         sd(adrenal_data$p4, na.rm = T),
         sd(adrenal_data$p5, na.rm = T),
         sd(adrenal_data$t, na.rm = T)),
  row.names = c('AE','AE2', 'B', 'DHC', 'DHT', 'E2', 'E3', 'P4', 'P5', 'T'))

liver_summary <- data.frame(
  mean = c(mean(liver_data$ae, na.rm = T),
           mean(liver_data$ae2, na.rm = T), 
           mean(liver_data$b, na.rm = T), 
           mean(liver_data$dhc, na.rm = T), 
           mean(liver_data$dht, na.rm = T), 
           mean(liver_data$e2, na.rm = T), 
           mean(liver_data$e3, na.rm = T), 
           mean(liver_data$p4, na.rm = T),
           mean(liver_data$p5, na.rm = T),
           mean(liver_data$t, na.rm = T)),
  n = c(length(liver_data$ae),
        length(liver_data$ae2), 
        length(liver_data$b), 
        length(liver_data$dhc), 
        length(liver_data$dht), 
        length(liver_data$e2), 
        length(liver_data$e3), 
        length(liver_data$p4),
        length(liver_data$p5),
        length(liver_data$t)),
  sd = c(sd(liver_data$ae, na.rm = T),
         sd(liver_data$ae2, na.rm = T), 
         sd(liver_data$b, na.rm = T), 
         sd(liver_data$dhc, na.rm = T), 
         sd(liver_data$dht, na.rm = T), 
         sd(liver_data$e2, na.rm = T), 
         sd(liver_data$e3, na.rm = T), 
         sd(liver_data$p4, na.rm = T),
         sd(liver_data$p5, na.rm = T),
         sd(liver_data$t, na.rm = T)),
  row.names = c('AE','AE2', 'B', 'DHC', 'DHT', 'E2', 'E3', 'P4', 'P5', 'T'))
#end of confusion - 4dec2022

#shapiro wilk test of normality  
raw_data %>% filter(steroid == "ae") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>% 
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "ae2") %>% 
  filter(matrix %in% c("adrenal")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "e2") %>% 
  filter(matrix %in% c("liver")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "e3") %>% 
  filter(matrix %in% c("adrenal")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "dht") %>% 
  filter(matrix %in% c("testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)
raw_data %>% filter(steroid == "p5") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  select(raw) %>% sapply(., shapiro.test)

#levene's test for homogeneity of variance 
#do you include matrix?? warning message about group coercion without it
raw_data %>% filter(steroid == "ae") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>% 
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "ae2") %>% 
  filter(matrix %in% c("adrenal")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "e2") %>% 
  filter(matrix %in% c("liver")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "e3") %>% 
  filter(matrix %in% c("adrenal")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "dht") %>% 
  filter(matrix %in% c("testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .) #dht only detectable in sti testes - no need for stats
raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)
raw_data %>% filter(steroid == "p5") %>% 
  filter(matrix %in% c("adrenal", "testes")) %>%
  leveneTest(raw ~ condition*matrix, data = .)


#summary tables
summarySE(data = raw_data, 
          measurevar = "raw", 
          groupvars = c("condition", "matrix", "steroid"))

raw_data %>% group_by(steroid, matrix, condition) %>% 
  dplyr::summarize(mean = mean(raw, na.rm = T),
                   sd = sd(raw, na.rm = T)) %>% 
  dplyr::arrange(matrix, steroid)

# androstenedione ####
#ae only detectable in both sti and control adrenal - so just a mann whitney on it
# add mann whitney stats to table
mw_ae <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_ae

ae <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% as.data.frame()

#boxplot
raw_data %>% filter(steroid == "ae") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "AE") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "ae") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "Androstenedione") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/ae_peripheral.pdf", width=6, height=4, dpi=100)

# alpha estradiol ####

#ae2 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_ae2 <- raw_data %>% 
  filter(steroid == "ae2") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = ., exact = FALSE)
mw_ae2

#boxplot
raw_data %>% filter(steroid == "ae2") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "aE2") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "ae2") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "alpha estradiol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .4), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/ae2_peripheral.pdf", width=6, height=4, dpi=100)


# estradiol ####
#e2 mann whitney liver - warning message: can't compute exact p-value with ties 
mw_e2 <- raw_data %>% 
  filter(steroid == "e2") %>% 
  filter(matrix == "liver") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e2

#boxplot
raw_data %>% filter(steroid == "e2") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "E2") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "e2") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "Estradiol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/e2_peripheral.pdf", width=6, height=4, dpi=100)


# estriol ####
#e3 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_e3 <- raw_data %>% 
  filter(steroid == "e3") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e3

#boxplot
raw_data %>% filter(steroid == "e3") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Estriol") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))


#barplot
raw_data %>% filter(steroid == "e3") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "estriol") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .2), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/e3_peripheral.pdf", width=6, height=4, dpi=100)


# corticosterone ####
#adrenals removed b/c cort concentrations exceed curve
res_aov_b <- raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_b) #effect of matrix, but not conditions


#pairwise posthoc
pwc_b_peripheral <-  raw_data %>% 
  filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_b_peripheral

#post-hoc with benjamini hochberg - main effect of sti on liver and pec, but not testes
one.way_b <- raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_b

#boxplot
raw_data %>% filter(steroid == "b") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "b") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "corticosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/b_peripheral.pdf", width=6, height=4, dpi=100)


# dehydrocorticosterone ####
#adrenals removed b/c dhc conc exceed curve
res_aov_dhc <- raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_dhc) #effect of matrix, but not condition

#pairwise posthoc
pwc_dhc_peripheral <-  raw_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_dhc_peripheral

#post-hoc with benjamini hochberg - main effect of sti on liver 
one.way_dhc <- raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_dhc

#boxplot
raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "dehydrocorticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))


#barplot
raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "dehydrocorticosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/dhc_peripheral.pdf", width=6, height=4, dpi=100)


# dihydrotestosterone ####
#dht only detectable in STI testes - no need for stats

#boxplot
raw_data %>% filter(steroid == "dht") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "dht") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "dihydrotestosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/dht_peripheral.pdf", width=6, height=4, dpi=100)


# testosterone ####
res_aov_t <- raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec", "adrenal", "testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_t)

#pairwise posthoc
pwc_t_peripheral <-  raw_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec","adrenal", "testes")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_t_peripheral

#post hoc looking at effect of sti on t in each matrix (not comparing across matrices)
one.way_t <- raw_data %>% filter(steroid == "t") %>% 
  filter(matrix %in% c("liver", "pec", "adrenal", "testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = log_con, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_t

#boxplot
raw_data %>% filter(steroid == "t") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "t") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "testosterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/t_peripheral.pdf", width=6, height=4, dpi=100)

# progesterone ####
#adrenals are removed b/c p4 concentration exceeds curve
res_aov_p4 <- raw_data %>% filter(steroid == "p4") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = log_con, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_p4)

#pairwise posthoc with benjamini hochberg correction - looks at differences across matrices
pwc_p4 <- raw_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_p4

#boxplot
raw_data %>% filter(steroid == "p4") %>%  filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "progesterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "progesterone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/p4_peripheral.pdf", width=6, height=4, dpi=100)


# pregnenolone ####

#run mann whitney b/c adrenals detectable in sti and control groups; testes only detectable in sti group
mw_p5 <- raw_data %>% 
  filter(steroid == "p5") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_p5


#boxplot
raw_data %>% filter(steroid == "p5") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "pregnenolone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))

#barplot
raw_data %>% filter(steroid == "p5") %>% 
  mutate(raw = coalesce(raw, 0)) %>%
  ggbarplot(x = "matrix", y = "raw",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "condition",
            palette = c("#bfbfbf", "#0E7C7B"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "pregnenolone") +
  theme(text = element_text(size = 18),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 120), expand = c(0, 0))
ggsave("~/Desktop/SG_juv_MS/r_figs/p5_peripheral.pdf", width=6, height=4, dpi=100)



