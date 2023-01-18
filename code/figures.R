library(tidyverse)
library(ggpubr) ; library(patchwork)
#library(kableExtra)


# Brain and Blood ####
# load data 
path <- c("SG_juv_MS/brain_blood/")
brain_data <- read_csv(paste0(path, "output_sofia_long_removed_wb.csv"))

path4figs <- c("SG_juv_MS/r_figs/")

## testosterone ####
t_fig_a <- 
  brain_data %>% filter(steroid == "t") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Testosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"#, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Testosterone",
            add = c("mean_se", "point"),
            width = .7,
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = NULL) +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .2), expand = c(0, 0))

t_fig_b <- 
  brain_data %>% filter(steroid == "t") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Testosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Testosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = NULL) +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0))

t_fig_a + t_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "t_brain.pdf"), 
       width = 12, height = 6, dpi = 100)

## corticosterone ####
b_fig_a <- 
  brain_data %>% filter(steroid == "b") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Corticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Corticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 40), expand = c(0, 0))

b_fig_b <- 
  brain_data %>% filter(steroid == "b") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Corticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Corticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0))

b_fig_a + b_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "b_brain.pdf"), 
       width = 12, height = 6, dpi = 100)

## dehydrocorticosterone #####
dhc_fig_a <-   
  brain_data %>% filter(steroid == "dhc") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>% 
  select(Matrix = matrix, Dehydrocorticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Dehydrocorticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0))

dhc_fig_b <-     
  brain_data %>% filter(steroid == "dhc") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>% 
  select(Matrix = matrix, Dehydrocorticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Dehydrocorticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .250001), expand = c(0, 0))

dhc_fig_a + dhc_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "dhc_brain.pdf"), 
       width = 12, height = 6, dpi = 100)  

## progesterone ####
p4_fig_a <- brain_data %>% 
  filter(steroid == "p4") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Progesterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Progesterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Concentration (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0))

p4_fig_b <- brain_data %>% 
  filter(steroid == "p4") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  select(Matrix = matrix, Progesterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "WB")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "WB", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Progesterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"), 
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = "") +
  theme(text = element_text(size = 15),
        legend.position = "left") +
  scale_y_continuous(limits = c(0, .300001), expand = c(0, 0))

p4_fig_a + p4_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "p4_brain.pdf"), 
       width = 12, height = 6, dpi = 100)

# all in one panel plot
t_fig_a + t_fig_b +
  b_fig_a + b_fig_b +
  dhc_fig_a + dhc_fig_b +
  p4_fig_a + p4_fig_b + 
  plot_layout(design = "ABBBBBB
                        CDDDDDD
                        EFFFFFF
                        GHHHHHH", 
              guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")


# Section jump ---------

# Behaviour ####

# load data 
path <- c("SG_juv_MS/behaviour/")
data <- read.csv(paste0(path, "sg_ms_behaviour.csv")) 
data %>% head()

# figures are made using raw data values (not log transformed)

# Collective latency boxplots ####
#note: adult controls - call latency missing for those that were noted to have called
#note: therefore inappropriate to include control adults in latency comparisons
data %>% dplyr::select(Condition = condition, Age = age, 
                       Call = call_lat, Song = song_lat,Flight = flight_lat,
                       Response = resp_lat) %>%
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult"),
         Condition = dplyr::recode(Condition, 
                                   sti = "STI", control = "Control")) %>%  
  gather(key, value, -Condition, - Age) %>% 
  mutate(ln_val = log1p(value)) %>% 
  ggplot(aes(x = Age, y = value, color = Condition)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .7, position=position_dodge(.75)) +
  labs(x = "", y = "") +
  scale_fill_brewer(palette="Dark2") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(4, "lines"), 
        strip.background = element_blank()) + 
  facet_wrap(~key, scales = "free", strip.position = "top", 
             labeller = labeller(key = c("Call"= "Call Latency", 
                                         "Flight"= "Flight Latency", 
                                         "Response"  = "Response Latency", 
                                         "Song" = "Song Latency")
             ))

## SONG LATENCY #####
#plot t-test for juv sti vs control - song latency boxplot 

data %>% filter(age=="juv") %>% 
  dplyr::select(Condition = condition, Song_Latency = song_lat) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>% 
  ggplot(aes(x = Condition, y = Song_Latency, fill = Condition)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  labs(x = "Condition", 
       y = "Song Latency (sec)", 
       subtitle = "Juvenile") +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t-test for juv sti vs control - song latency bar graph 

#following code uses ggplot with geom_bar 
#you have to manually calculate mean and SEM; can't have individual data points
data %>% filter(age=="juv") %>% 
  dplyr::select(Condition = condition, Song_Latency = song_lat) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>%
  group_by(Condition) %>%
  summarize(aver = mean(Song_Latency), 
            n = n(Condition), se = sd(Song_Latency)/sqrt(n)) %>% 
  ggplot(aes(x = Condition, y = aver, fill = Condition)) +
  geom_bar(position = position_dodge(), stat = "identity", width = .9)  +
  geom_errorbar(aes(ymin=aver-se, ymax=aver+se),
                width = .1, lwd = .8,                     # Width of the error bars
                position=position_dodge(.9)) + 
  labs(x = "Condition", 
       y = "Song Latency (seconds)",
       subtitle = "Juvenile") +
  ggthemes::theme_tufte(base_size = 15) +
  scale_fill_grey(start = .7, end = .4) + 
  theme(panel.spacing = unit(1, "lines"), 
        strip.background = element_blank(),  
        legend.position = "none",
        axis.line = element_line(size = .5, colour = "black", linetype=1)) +
  scale_y_continuous(limits = c(0, 650), expand = c(0, 0))

#following code uses ggpubr with ggbarplot - you don't have to manually calculate mean and SEM, but I can't figure out how to adjust size of individual values without adjusting error bars, jitter intensity (either use "jitter" or "dotplot" in add command)

#song_juv <- 
data %>% filter(age=="juv") %>% 
  select(Condition = condition, Song_Latency = song_lat) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>%
  ggbarplot(x = "Condition", y = "Song_Latency", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#4E606A"),
            position = position_dodge(1)) +
  theme(legend.position = "none") + 
  labs(x = "Condition", 
       y = "Song Latency (seconds)", 
       subtitle = "Juvenile") +
  scale_y_continuous(limits = c(0, 650), expand = c(0, 0))

# plot t-test for STI adults vs juveniles - song latency boxplot
data %>% filter(condition=="sti") %>% 
  select(Age = age, Song_Latency = song_lat) %>% 
  mutate(Age = recode(Age, juv = "Juveniles", adult = "Adult")) %>% 
  ggplot(aes(x = Age, y = Song_Latency, fill = Age)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  labs(x = "Age", 
       y = "Song Latency (sec)", 
       subtitle = "STI") +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

## plot t-test for STI adults vs juveniles - song latency bar graph 
data %>% filter(condition=="sti") %>% 
  select(Age = age, Song_Latency = song_lat) %>% 
  mutate(Age = recode(Age, juv = "Juveniles", adult = "Adult")) %>% 
  ggbarplot(x = "Age", y = "Song_Latency", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Age",
            palette = c("#bfbfbf", "#4E606A"),
            position = position_dodge(1)) +
  theme(legend.position = "none") + 
  labs(x = "Age", 
       y = "Song Latency (seconds)", 
       subtitle = "STI") +
  scale_y_continuous(limits = c(0, 650), expand = c(0, 0))


## RESPONSE LATENCY #####
#plot t test for juv sti vs control - response latency boxplot 
data %>% filter(age=="juv") %>% 
  select(Condition = condition, Response_Latency = resp_lat) %>% 
  mutate(Condition = recode(Condition, sti = "STI", control = "Control")) %>% 
  ggplot(aes(x = Condition, y = Response_Latency, fill = Condition)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Condition", 
       y = "Response Latency (sec)", 
       subtitle = "Juvenile") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for juv sti vs control - response latency bar graph
juv_resp_lat <- 
  data %>% filter(age == "juv") %>% 
  select(Condition = condition, Response_Latency = resp_lat) %>% 
  mutate(Condition = dplyr::recode(Condition, control = "Control", sti = "STI")) %>% 
  ggbarplot(x = "Condition", y = "Response_Latency", 
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1),
            color = "black", fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Response Latency (sec)"
       #,subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 606.0000001), expand = c(0, 0))

#plot t test for STI adults vs juveniles - response latency boxplot
data %>% filter(condition=="sti") %>% 
  select(Age = age, Response_Latency = resp_lat) %>% 
  mutate(Age = recode(Age, juv = "Juveniles", adult = "Adult")) %>% 
  ggplot(aes(x = Age, y = Response_Latency, fill = Age)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Age", 
       y = "Response Latency (sec)", 
       subtitle = "STI") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for STI adults vs juveniles - response latency bar graph
sti_resp_lat <- 
  data %>% filter(condition == "sti") %>% 
  select(Age = age, Response_Latency = resp_lat) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggbarplot(x = "Age", y = "Response_Latency", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Age",
            palette = c("#505053", "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Response Latency (sec)" 
       #, subtitle = "STI"
  ) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0))

## NUMBER OF SONGS #####
#number of songs boxplot ANOVA
data %>% select(Condition = condition, Age = age, 
                Songs = songs) %>%
  mutate(Age = recode(Age, juv = "Juvenile", adult = "Adult"),
         Condition = recode(Condition, sti = "STI", control = "Control")) %>%  
  gather(key, value, -Condition, - Age) %>% 
  ggplot(aes(x = Age, y = value, color = Condition)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  labs(x = "", y = "") +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Condition", 
       y = "Number of Songs") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank()) + 
  facet_wrap(~key, scales = "free", strip.position = "top", 
             labeller = labeller(key = c("Songs"= "Number of Songs", 
                                         "Flights"= "Number of Flights", 
                                         "Time" = "Time in 5M")
             ))


#plot t test for juv sti vs control - number of songs boxplot
data %>% filter(age=="juv") %>% 
  select(Condition = condition, Songs = songs) %>% 
  mutate(Condition = recode(Condition, sti = "STI", control = "Control")) %>% 
  ggplot(aes(x = Condition, y = Songs, fill = Condition)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Condition", 
       y = "Number of Songs", 
       subtitle = "Juvenile") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for juv sti vs control - number of songs bar graph
juv_n_song <- 
  data %>% filter(age == "juv") %>% 
  select(Condition = condition, Songs = songs) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>% 
  ggbarplot(x = "Condition", y = "Songs",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Number of Songs" 
       #,subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

#plot t test for STI adults vs juveniles - number of songs boxplot
data %>% filter(condition=="sti") %>% 
  select(Age = age, Songs = songs) %>% 
  mutate(Age = recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggplot(aes(x = Age, y = Songs, fill = Age)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Age", 
       y = "Number of Songs", 
       subtitle = "STI") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for STI adults vs juveniles - number of songs bar graph
sti_n_song <- 
  data %>% filter(condition == "sti") %>% 
  select(Age = age, Songs = songs) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>%
  ggbarplot(x = "Age", y = "Songs", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Age",
            palette = c("#505053", "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Number of Songs", 
       #subtitle = ""
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))

## NUMBER OF FLIGHTS ####
#number of flights boxplot ANOVA
data %>% select(Condition = condition, Age = age, 
                Flights = flights) %>%
  mutate(Age = recode(Age, juv = "Juvenile", adult = "Adult"),
         Condition = recode(Condition, sti = "STI", control = "Control")) %>%  
  gather(key, value, -Condition, - Age) %>% 
  ggplot(aes(x = Age, y = value, color = Condition)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  labs(x = "", y = "") +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Condition", 
       y = "Number of Flights") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank()) + 
  facet_wrap(~key, scales = "free", strip.position = "top", 
             labeller = labeller(key = c("Songs"= "Number of Songs", 
                                         "Flights"= "Number of Flights", 
                                         "Time" = "Time in 5M")
             ))

#plot t test for juv sti vs control - number of flights boxplot
data %>% filter(age=="juv") %>% 
  select(Condition = condition, Flights = flights) %>% 
  mutate(Condition = recode(Condition, sti = "STI", control = "Control")) %>% 
  ggplot(aes(x = Condition, y = Flights, fill = Condition)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Condition", 
       y = "Number of Flights", 
       subtitle = "Juvenile") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for juv sti vs control - number of flights bargraph
juv_fl <- 
  data %>% filter(age == "juv") %>% 
  select(Condition = condition, Flights = flights) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>% 
  ggbarplot(x = "Condition", y = "Flights",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Number of Flights", 
       #subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0))

#plot t test for STI adults vs juveniles - number of flights boxplot
data %>% filter(condition=="sti") %>% 
  select(Age = age, Flights = flights) %>% 
  mutate(Age = recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggplot(aes(x = Age, y = Flights, fill = Age)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Age", 
       y = "Number of Flights", 
       subtitle = "STI") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for STI adults vs juveniles - number of flights bargraph
sti_fl <- 
  data %>% filter(condition == "sti") %>% 
  select(Age = age, Flights = flights) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggbarplot(x = "Age", y = "Flights", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Age",
            palette = c("#505053", "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Number of Flights", 
       #subtitle = "STI"
  ) +
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0))


## TIME IN 5 METERS ####
#time spent within 5M boxplot ANOVA
data %>% select(Condition = condition, Age = age, 
                Time = time) %>%
  mutate(Age = recode(Age, juv = "Juvenile", adult = "Adult"),
         Condition = recode(Condition, sti = "STI", control = "Control")) %>%  
  gather(key, value, -Condition, - Age) %>% 
  ggplot(aes(x = Age, y = value, color = Condition)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  labs(x = "", y = "") +
  scale_fill_brewer(palette="Dark2") +
  labs(x = "Condition", 
       y = "Time in 5M (sec)") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank()) + 
  facet_wrap(~key, scales = "free", strip.position = "top", 
             labeller = labeller(key = c("Songs"= "Number of Songs", 
                                         "Flights"= "Number of Flights", 
                                         "Time" = "Time in 5M")
             ))

#plot t test for juv sti vs control - time in 5M boxplot
data %>% filter(age=="juv") %>% 
  select(Condition = condition, Time_in_5M = time) %>% 
  mutate(Condition = recode(Condition, sti = "STI", control = "Control")) %>% 
  ggplot(aes(x = Condition, y = Time_in_5M, fill = Condition)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Condition", 
       y = "Time Spent in 5M (seconds)", 
       subtitle = "Juvenile") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for juv sti vs control - time in 5M bar graph
juv_t <- 
  data %>% filter(age=="juv") %>% 
  select(Condition = condition, Time_in_5M = time) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>%  
  ggbarplot(x = "Condition", y = "Time_in_5M",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Condition",
            palette = c("#bfbfbf", "#505053"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Time Spent in 5m (sec)", 
       #subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 605), expand = c(0, 0))

#plot t test for STI adults vs juveniles - time in 5M boxplot
data %>% filter(condition=="sti") %>% 
  select(Age = age, Time_in_5M = time) %>% 
  mutate(Age = recode(Age, juv = "Juveniles", adult = "Adult")) %>% 
  ggplot(aes(x = Age, y = Time_in_5M, fill = Age)) +
  geom_boxplot(color = "black", alpha = 10) +
  geom_dotplot(binaxis='y', stackdir='center', 
               dotsize = .5, position=position_dodge(.75)) +
  scale_fill_manual(values = c("#4863A0", "#8C001A")) +
  labs(x = "Age", 
       y = "Time in 5M (seconds)", 
       subtitle = "STI") +
  theme_classic(base_size = 17) +
  theme(panel.spacing = unit(2, "lines"), 
        strip.background = element_blank(), 
        legend.position = "none")

#plot t test for STI adults vs juveniles - time in 5M bar graph
sti_t <- 
  data %>% filter(condition=="sti") %>% 
  select(Age = age, Time_in_5M = time) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>%   
  ggbarplot(x = "Age", y = "Time_in_5M", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .5, alpha = 1), 
            fill = "Age",
            palette = c("#505053", "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 18)) + 
  labs(x = "", 
       y = "Time Spent in 5m (sec)", 
       #subtitle = "STI"
  ) +
  scale_y_continuous(limits = c(0, 608.0000001), expand = c(0, 0))

## panel plots ####
juv_resp_lat + juv_t + juv_n_song + juv_fl + 
  plot_layout(design = "AB
  						CD", 
              guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")
ggsave(paste0(path4figs, "behav_juv.pdf"), 
       width = 12, height = 12, dpi = 100)


sti_resp_lat + sti_t + sti_n_song + sti_fl +
  plot_layout(design = "AB
  						CD", 
              guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")
ggsave(paste0(path4figs, "behav_sti.pdf"), 
       width = 12, height = 12, dpi = 100)
