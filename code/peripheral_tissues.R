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

# normality test shapiro wilk  ####
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

# levene's test for homogeneity of variance ####
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

# summary tables #####
summarySE(data = raw_data, 
          measurevar = "raw", 
          groupvars = c("condition", "matrix", "steroid"))

raw_data %>% 
  group_by(steroid, matrix, condition) %>% 
  dplyr::summarize(mean = mean(raw, na.rm = T),
                   sd = sd(raw, na.rm = T)) %>% 
  dplyr::arrange(matrix, steroid)

raw_data %>% 
  filter(!is.na(log_con)) %>% 
  group_by(steroid, matrix) %>% 
  tally() %>% view()


## androstenedione ####
#ae only detectable in both sti and control adrenal - so just a mann whitney on it

ae.out <- raw_data %>% 
  filter(steroid == "ae") %>% 
  lmer(log_con ~ condition + (1|subject), data = .)
anova(ae.out) 

#boxplot
raw_data %>% filter(steroid == "ae") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "AE") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype = 1)) +
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
ggsave("~/Desktop/SG_juv_MS/r_figs/ae_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## alpha estradiol ####

ae2.out <- raw_data %>% 
  filter(steroid == "ae2") %>% 
  lmer(log_con ~ condition + (1|subject), data = .)
anova(ae2.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/ae2_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## estradiol ####
e2.out <- raw_data %>% 
  filter(steroid == "e2") %>% 
  lmer(log_con ~ condition + (1|subject), data = .)
anova(e2.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/e2_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## estriol ####

e3.out <- raw_data %>% 
  filter(steroid == "e3") %>% 
  lmer(log_con ~ condition + (1|subject), data = .)
anova(e3.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/e3_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## corticosterone ####

b.out <- raw_data %>% 
  filter(steroid == "b") %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(b.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/b_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## dehydrocorticosterone ####

dhc.out <- raw_data %>% 
  filter(steroid == "dhc") %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(dhc.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/dhc_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## dihydrotestosterone ####
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
ggsave("~/Desktop/SG_juv_MS/r_figs/dht_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)


## testosterone ####

tt.out <- raw_data %>% 
  filter(steroid == "t") %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(tt.out) 

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
ggsave("~/Desktop/SG_juv_MS/r_figs/t_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

## progesterone ####

p4.out <- raw_data %>% 
  filter(steroid == "p4") %>% 
  lmer(log_con ~ condition*matrix + (1|subject), data = .)
anova(p4.out) 

#boxplot
raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
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
ggsave("~/Desktop/SG_juv_MS/r_figs/p4_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)


## pregnenolone ####

p5.out <- raw_data %>% 
  filter(steroid == "p5") %>% 
  lmer(log_con ~ condition + (1|subject), data = .)
anova(p5.out) 


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
ggsave("~/Desktop/SG_juv_MS/r_figs/p5_peripheral.pdf", 
       width = 6, height = 4, dpi = 100)

# #####
sjPlot::tab_model(ae.out, ae2.out, b.out, dhc.out, tt.out, p4.out, p5.out, 
                  p.val = "kr", collapse.ci = TRUE, 
                  show.reflvl = F, prefix.labels = "varname", 
                  #, p.style = "stars"
                  title = "Periferal tissues"
                  , dv.labels = c("ae", "ae2",
                                  "b", "dhc",
                                  "t", "p4", "p5")
                  , file = paste0(path4figs, "lmm_sum_perif.html"))
webshot::webshot(paste0(path4figs, "lmm_sum_perif.html"), 
                 paste0(path4figs, "lmm_sum_perif.png"), zoom = 2)
# LMM - diagnostics #### 
library(DHARMa)
md_ls <- list(ae.out, ae2.out, b.out, dhc.out, tt.out, p4.out, p5.out)
names(md_ls) <- c("ae.out", "ae2.out", "b.out", "dhc.out", "tt.out", "p4.out", "p5.out")

path4figs <- c("SG_juv_MS/r_figs/")

for(i in 1:length(md_ls)){
  model2diag <- md_ls[[i]]
  pdf(paste0(path4figs, "udoz_perif_", names(md_ls)[i],".pdf"), 
      onefile = T, width = 7, height = 7, compress = T)
  test_res <- testResiduals(model2diag) ## uniformity, dispersion, and outliers test
  testZeroInflation(model2diag)
  dev.off()
}


# Other models ####

## androstenedione ####
# add mann whitney stats to table
mw_ae <- raw_data %>% 
  filter(steroid == "ae") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_ae

## alpha estradiol ####

#ae2 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_ae2 <- raw_data %>% 
  filter(steroid == "ae2") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = ., exact = FALSE)
mw_ae2

# #estradiol ####
#e2 mann whitney liver - warning message: can't compute exact p-value with ties 
mw_e2 <- raw_data %>% 
  filter(steroid == "e2") %>% 
  filter(matrix == "liver") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e2

## estriol ####
#e3 mann whitney adrenals - warning message: can't compute exact p-value with ties 
mw_e3 <- raw_data %>% 
  filter(steroid == "e3") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_e3

# #corticosterone ####
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

## dehydrocorticosterone ####
#adrenals removed b/c dhc conc exceed curve
res_aov_dhc <- raw_data %>% filter(steroid == "dhc") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
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

## testosterone ####
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

## progesterone ####
#adrenals are removed b/c p4 concentration exceeds curve
res_aov_p4 <- raw_data %>% filter(steroid == "p4") %>% 
  filter(matrix %in% c("liver", "pec","testes")) %>% 
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

## pregnenolone ####

#run mann whitney b/c adrenals detectable in sti and control groups; testes only detectable in sti group
mw_p5 <- raw_data %>% 
  filter(steroid == "p5") %>% 
  filter(matrix == "adrenal") %>% 
  wilcox.test(raw ~ condition, data = .)
mw_p5
