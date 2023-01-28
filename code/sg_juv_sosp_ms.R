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

#do NOT run stats with NAs as zero !!

## testosterone ####
# anova is throwing error b/c AH, VMH, VTA, BNST, and Cb is detectable ONLY 
# in STI or control group - should I exclude from analyses then?  
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

# I know I shouldn't do these b/c no significant effects from ANOVA table- just to set up code
#pairwise posthoc: sti vs control grouped by matrix
pwc_t <-  brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%   
  group_by(matrix) %>%
  pairwise_t_test(log_con ~ condition, p.adjust.method = "BH")
pwc_t

#lmm: also can't handle that some matrices have NA for an entire condition; lmm auto-drops these columns
#here, interaction is significant (but sketchy); F value not that different btw ANOVA and LMM - maybe correction is changing p value
t.brain.out <- brain_data %>% 
  filter(steroid == "t") %>% 
  filter(matrix %in% c("ah", "bnst", "cb", "cg", 
                       "ls", "ncm", "poa", "tna", "vmh", "vta")) %>% 
  lmer(log_con ~ matrix*condition + (1|subject), data = .)
anova(t.brain.out)

# lmm post-hoc pairwise comparison
t.brain.emm <- 
  emmeans::emmeans(t.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")

t.brain.emm$`pairwise differences of matrix, condition` %>% 
  data.frame() %>% 
  filter(p.value <= .1) %>% 
  mutate(p.value = round(p.value, 3)) %>% 
  separate(X1, c("pair_1", "pair_2"), sep = " - ") %>% 
  separate(pair_1, c("m_1", "con_1"), sep = " ") %>% 
  separate(pair_2, c("m_2", "con_2"), sep = " ") %>%
  filter(!con_1==con_2) %>% 
  kbl(caption = "Brain - b") %>%
  kable_paper(bootstrap_options = "striped", 
              full_width = F, html_font = 12) %>% 
  save_kable(file = paste0(path4figs, "lmm_posthoc_brain_t.html"), 
             self_contained = T)

# mann whitney on blood testosterone
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
b.brain.out <- brain_data %>% 
    filter(steroid == "b") %>% 
    filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls",
                         "ncm", "poa", "tna", "vmh", "vta")) %>% 
    lmer(log_con ~ matrix*condition + (1|subject), data = .)
anova(b.brain.out) 

# lmm post-hoc pairwise comparison
b.brain.emm <- 
  emmeans::emmeans(b.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")

b.brain.emm$`pairwise differences of matrix, condition` %>% 
  data.frame() %>% 
  filter(p.value <= .1) %>% 
  mutate(p.value = round(p.value, 3)) %>% 
  separate(X1, c("pair_1", "pair_2"), sep = " - ") %>% 
  separate(pair_1, c("m_1", "con_1"), sep = " ") %>% 
  separate(pair_2, c("m_2", "con_2"), sep = " ") %>%
  filter(!con_1==con_2) %>% 
  kbl(caption = "Brain - b") %>%
  kable_paper(bootstrap_options = "striped", 
              full_width = F, html_font = 12) %>% 
  save_kable(file = paste0(path4figs, "lmm_posthoc_brain_b.html"), 
             self_contained = T)
  
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
dhc.brain.out <- brain_data %>% 
  filter(steroid == "dhc") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", 
                       "ncm", "poa", "tna", "vmh", "vta")) %>% 
lmer(log_con ~ condition * matrix + (1|subject), data = .)
dhc.brain.out %>% anova()

# lmm post-hoc pairwise comparison
dhc.brain.emm <- 
  emmeans::emmeans(dhc.brain.out, 
                   list(pairwise ~  matrix * condition),
                   adjust = "bh", lmer.df = "satterthwaite")

dhc.brain.emm$`pairwise differences of matrix, condition` %>% 
  data.frame() %>% 
  filter(p.value <= .1) %>% 
  mutate(p.value = round(p.value, 3)) %>% 
  separate(X1, c("pair_1", "pair_2"), sep = " - ") %>% 
  separate(pair_1, c("m_1", "con_1"), sep = " ") %>% 
  separate(pair_2, c("m_2", "con_2"), sep = " ") %>%
  filter(!con_1==con_2) %>% 
  kbl(caption = "Brain - dhc") %>%
  kable_paper(bootstrap_options = "striped", 
              full_width = F, html_font = 12) %>% 
  save_kable(file = paste0(path4figs, "lmm_posthoc_brain_dhc.html"), 
             self_contained = T)

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
  
# lmm - again, I dropped condition as a fixed effect b/c p4 only detectable in STI group
p4.brain.out <- brain_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("ah" , "bnst", "cb", "cg", "ls", 
                       "ncm", "poa", "tna", "vmh", "vta")) %>% 
lmer(log_con ~ matrix + (1|subject), data = .)
anova(p4.brain.out)

# lmm post-hoc pairwise comparison
p4.brain.emm <- 
  emmeans::emmeans(p4.brain.out, 
                   list(pairwise ~  matrix),
                   adjust = "bh", lmer.df = "satterthwaite")

p4.brain.emm$`pairwise differences of matrix` %>% 
  data.frame() %>% 
  filter(p.value <= .1) %>% 
  kbl(caption = "Brain - p4") %>%
  kable_paper(bootstrap_options = "striped", 
              full_width = F, html_font = 12) %>% 
  save_kable(file = paste0(path4figs, "lmm_posthoc_brain_p4.html"), 
             self_contained = T)

# mann whitney on blood progesterone
mw_p4_circ <- brain_data %>% 
  filter(steroid == "p4") %>% 
  filter(matrix %in% c("wb")) %>% 
  wilcox.test(concentration ~ condition, data = .)
mw_p4_circ

# ####
sjPlot::tab_model(t.brain.out, b.brain.out, dhc.brain.out, p4.brain.out, 
                  p.val = "kr", collapse.ci = TRUE, 
                  show.reflvl = F, prefix.labels = "varname", 
                  #, p.style = "stars"
                  title = "Brain"
                  , dv.labels = c("Testosterone", 
                                  "Corticosterone",
                                  "Dehydrocorticosterone",
                                  "Progesterone")
                  , file = paste0(path4figs, "lmm_sum_brain.html"))
  
webshot::webshot(paste0(path4figs, "lmm_sum_brain.html"), 
                 paste0(path4figs, "lmm_sum_brain.png"), zoom = 2)

# LMM - diagnostics #### 
library(DHARMa)
md_ls <- list(t.brain.out, b.brain.out, dhc.brain.out, p4.brain.out)
names(md_ls) <- c("brain.t.out", "brain.b.out", "brain.dhc.out", "brain.p4.out")

path4figs <- c("SG_juv_MS/r_figs/")

for(i in 1:length(md_ls)){
  model2diag <- md_ls[[i]]
  pdf(paste0(path4figs, "udoz_", names(md_ls)[i], ".pdf"), 
      onefile = T, width = 7, height = 7, compress = T)
  test_res <- testResiduals(model2diag) ## uniformity, dispersion, and outliers test
  testZeroInflation(model2diag)
  dev.off()
 }


