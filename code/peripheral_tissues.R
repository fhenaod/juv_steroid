#peripheral tissues in juvenile song sparrows
#7 November 2022

#the raw column is ng/g (the other columns are incorrect) 
#imputations were made on ng/g values, NOT pg values

setwd("~/Desktop/SG_juv_MS/peripheral_tissue/")
raw_data <- read.csv("output_data.csv") %>% select(-tissue_weight, -conc)
library(tidyverse)
library(rstatix)

raw_data %>% ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1)
raw_data %>% filter(steroid == "ae") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "ae2") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "e2") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "e3") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "b") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "dhc") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "dht") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "t") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "p4") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()
raw_data %>% filter(steroid == "p5") %>% ggplot(aes(x = condition, y = raw)) +
  geom_boxplot()


# androstenedione anova and fig by matrix ####
res_aov_ae <- raw_data %>% filter(steroid == "ae") %>% 
  anova_test(data = ., dv = raw, wid = subject, 
             between = condition)
get_anova_table(res_aov_ae)

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

# alpha estradiol anova and fig by matrix ####
res_aov_ae2 <- raw_data %>% filter(steroid == "ae2") %>% 
  anova_test(data = ., dv = raw, wid = subject,
             between = condition)
get_anova_table(res_aov_ae2)

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

# estradiol anova and fig by matrix ####
res_aov_e2 <- raw_data %>% filter(steroid == "e2") %>% 
  anova_test(data = ., dv = raw, wid = subject, between = condition)
get_anova_table(res_aov_e2)

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

# estriol anova and fig by matrix ####
res_aov_e3 <- raw_data %>% filter(steroid == "e3") %>% 
  anova_test(data = ., dv = raw, wid = subject,
             between = condition)
get_anova_table(res_aov_e3)

raw_data %>% filter(steroid == "e3") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "E3") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, .5), expand = c(0, 0))

# corticosterone anova and fig by matrix ####
#adrenals removed b/c cort concentrations exceed curve
res_aov_b <- raw_data %>% filter(steroid == "b") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = raw, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_b) #effect of matrix, but not conditions

#post-hoc with benjamini hochberg - main effect of sti on liver and pec, but not testes
one.way_b <- raw_data %>% filter(steroid == "b") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  group_by(matrix) %>%
  anova_test(dv = raw, wid = subject, between = condition) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
one.way_b

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

# dehydrocorticosterone anova and fig by matrix ####
#adrenals removed b/c dhc conc exceed curve
res_aov_dhc <- raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = raw, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_dhc) #effect of matrix, but not condition

raw_data %>% filter(steroid == "dhc") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))

# dihydrotestosterone anova and fig by matrix ####
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

# testosterone anova and fig by matrix ####
res_aov_t <- raw_data %>% filter(steroid == "t") %>% 
  anova_test(data = ., dv = raw, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_t)

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

# progesterone anova and fig by matrix ####
#adrenals are removed b/c p4 concentration exceeds curve
res_aov_p4 <- raw_data %>% filter(steroid == "p4") %>% filter(matrix %in% c("liver", "pec","testes")) %>% 
  anova_test(data = ., dv = raw, wid = subject,
             within = matrix, between = condition)
get_anova_table(res_aov_p4)

raw_data %>% filter(steroid == "p4") %>%  filter(matrix %in% c("liver", "pec","testes")) %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0))

# pregnenolone anova and fig by matrix ####
res_aov_p5 <- raw_data %>% filter(steroid == "p5") %>% 
  anova_test(data = ., dv = raw, wid = subject,
             between = condition)
get_anova_table(res_aov_p5)

raw_data %>% filter(steroid == "p5") %>% 
  ggplot(aes(x = matrix, y = raw, fill = condition)) +
  geom_boxplot(color = "black", alpha = 1) + theme_classic() + 
  labs(x = "Matrix", 
       y = "ng/g or ng/ml",
       subtitle = "Corticosterone") +
  scale_fill_grey(start = .9, end = .4) +
  theme(panel.spacing = unit(1, "lines"), 
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))


