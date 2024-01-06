library(tidyverse)
library(ggpubr) ; library(patchwork)
library(wesanderson)
library(ggbreak)
#library(kableExtra)


# Brain and Blood ####
# load data 
brain_data <- read.csv("rev2/sg_ms_brain.csv")
path4figs <- c("rev2/")

## testosterone ####
t_fig_a <- 
  brain_data %>% filter(steroid == "t") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Testosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"#, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Testosterone",
            add = c("mean_se", "point"),
            width = .7,
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Testosterone (ng/g)", 
       subtitle = NULL) +
  theme(text = element_text(size = 18), 
        legend.position = "none", legend.title = element_blank()) +
  scale_y_continuous(limits = c(0, .16), expand = c(0, 0))

t_fig_b <- 
  brain_data %>% filter(steroid == "t") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Testosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Testosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = NULL) +
  theme(legend.position = "none", text = element_text(size = 18), 
        #legend.title = element_blank(), 
        #legend.position = c(0, 0.2), 
        #legend.justification = c(-0.15, 0),
        #legend.direction = "vertical"
        ) +
  scale_y_continuous(limits = c(0, 3.25), 
                     expand = c(0, 0), breaks = c(0, .25, .5, .75, 1, 3.25)) +
  geom_text(x = 1.82, y = 0.05, label = "nd", size = 6) +
  geom_text(x = 4.18, y = 0.05, label = "nd", size = 6) +
  geom_text(x = 4.82, y = 0.05, label = "nd", size = 6) +
  geom_text(x = 6.18, y = 0.05, label = "nd", size = 6) +
  geom_text(x = 9.82, y = 0.05, label = "nd", size = 6) +
    scale_y_cut(breaks = c(1.1), which = c(1, 2), 
                scales = c(0.1, 1), space = c(1.25)) 

t_fig_a + t_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'keep') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "t_brain.pdf"), 
       width = 16, height = 8, dpi = 300)
ggsave(paste0(path4figs, "t_brain.svg"), 
       width = 16, height = 8, dpi = 300)

## corticosterone ####
b_fig_a <- 
  brain_data %>% filter(steroid == "b") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Corticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Corticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Corticosterone (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
  geom_signif(#comparisons = list(c("Control", "STI")), 
    annotations = "****",
    y_position = 15.6, xmin = .82, xmax = 1.18, 
    tip_length = 10, vjust = .5, size = 1, textsize = 8)

b_fig_b <- 
  brain_data %>% filter(steroid == "b") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Corticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Corticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(), 
        legend.position = c(0, 0.85), 
        legend.justification = c(-0.15, 0),
        legend.direction = "vertical") +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
  geom_signif(annotations = c("****", "****", "***", "***", "****", "****", "****", "***", "*", "****"), 
              y_position = c(5.5, 5.5, 5.5, 2.6, 5.5, 5.8, 5.5, 6.6, 5.5, 6.9), 
              xmin = seq(from = .82, to = 9.82, by = 1), 
              xmax = seq(from = 1.18, to = 10.18, by = 1), 
              tip_length = 10, vjust = .5, size = 1, textsize = 8)

b_fig_a + b_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'keep') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "b_brain.pdf"), 
       width = 16, height = 8, dpi = 300)
ggsave(paste0(path4figs, "b_brain.svg"), 
       width = 16, height = 8, dpi = 300)

## dehydrocorticosterone #####
dhc_fig_a <-   
  brain_data %>% filter(steroid == "dhc") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>% 
  dplyr::select(Matrix = matrix, Dehydrocorticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Dehydrocorticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Dehydrocorticosterone (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  geom_signif(#comparisons = list(c("Control", "STI")), 
              annotations = "****",
              y_position = 3.6, xmin = .82, xmax = 1.18, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

dhc_fig_b <-     
  brain_data %>% filter(steroid == "dhc") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>% 
  dplyr::select(Matrix = matrix, Dehydrocorticosterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Dehydrocorticosterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(), 
        legend.position = c(0, 0.85), 
        legend.justification = c(-0.15, 0),
        legend.direction = "vertical") +
  scale_y_continuous(limits = c(0, .250001), expand = c(0, 0)) + 
  geom_signif(annotations = c("#", "#"), 
              y_position = c(0.14, 0.14), 
              xmin = c(0.82, 3.82), xmax = c(1.18, 4.18), 
              tip_length = 10, vjust = -.2, size = 1, textsize = 8) + 
  geom_text(x = 9.82, y = 0.004, label = "nd", size = 6)
  
dhc_fig_a + dhc_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'keep'
              ) +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "dhc_brain.pdf"), 
       width = 16, height = 8, dpi = 300)  

## progesterone ####
p4_fig_a <- brain_data %>% 
  filter(steroid == "p4") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Progesterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("wb"
                       #, "plasma"
  )) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Progesterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = "Progesterone (ng/g)", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
  geom_signif(#comparisons = list(c("Control", "STI")), 
    annotations = "#",
    y_position = 2.65, xmin = .82, xmax = 1.18, 
    tip_length = 10, vjust = -.2, size = 1, textsize = 8)

p4_fig_b <- brain_data %>% 
  filter(steroid == "p4") %>% 
  mutate(concentration = coalesce(concentration, 0)) %>%
  dplyr::select(Matrix = matrix, Progesterone = concentration, Condition = condition) %>% 
  filter(Matrix %in% c("ah" , "bnst", "cb", "cg", "ls", "ncm", "poa", "tna", "vmh", "vta")) %>%  
  mutate(Matrix = dplyr::recode(Matrix, ah = "AH", bnst = "BnST", cb = "Cb", cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", plasma = "Plasma", poa = "POA", tna = "TnA", vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", con = "Control")) %>% 
  mutate(Matrix = fct_relevel(Matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>%
  ggbarplot(x = "Matrix", y = "Progesterone",
            add = c("mean_se", "point"), 
            add.params = list(size = .8, alpha = 1), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.7)) +
  labs(x = NULL, 
       y = " ", 
       subtitle = "") +
  theme(text = element_text(size = 18), legend.title = element_blank(), 
        legend.position = c(0, 0.85), 
        legend.justification = c(-0.15, 0),
        legend.direction = "vertical") +
  scale_y_continuous(limits = c(0, .300001), expand = c(0, 0)) +
  geom_text(x = 0.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 1.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 2.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 3.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 4.18, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 4.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 5.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 6.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 7.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 8.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 9.82, y = 0.005, label = "nd", size = 6) +
  geom_text(x = 10.18,y = 0.005, label = "nd", size = 6)

p4_fig_a + p4_fig_b + 
  plot_layout(design = "ABBBBBB", guides = 'keep') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")

ggsave(paste0(path4figs, "p4_brain.pdf"), 
       width = 16, height = 8, dpi = 300)

# all in one panel plot
t_fig_a + t_fig_b +
  b_fig_a + b_fig_b +
  dhc_fig_a + dhc_fig_b +
  p4_fig_a + p4_fig_b + 
  plot_layout(design = "ABBBBBB
                        CDDDDDD
                        EFFFFFF
                        GHHHHHH", 
              guides = 'keep') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")


# Section jump ---------

# Behaviour ####

# load data 
path <- c("data/")
data <- read.csv(paste0(path, "sg_ms_behaviour.csv")) 
data %>% head()

# figures are made using raw data values (not log transformed)

## RESPONSE LATENCY #####
#plot t test for juv sti vs control - response latency bar graph
library(wesanderson)
juv_resp_lat <- 
  data %>% filter(age == "juv") %>% 
  dplyr::select(Condition = condition, Response_Latency = resp_lat) %>% 
  mutate(Condition = dplyr::recode(Condition, control = "Control", sti = "STI")) %>% 
  ggbarplot(x = "Condition", y = "Response_Latency", 
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"),
            color = "black", fill = "Condition",
           #palette = c("#bfbfbf", "#505053"), # gray scale
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(.1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Response Latency (sec)"
       #,subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 660), expand = c(0, 0)) +
  geom_signif(comparisons = list(c("Control", "STI")), annotations = "****",
              y_position = 600, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

#plot t test for STI adults vs juveniles - response latency bar graph
sti_resp_lat <- 
  data %>% filter(condition == "sti") %>% 
  dplyr::select(Age = age, Response_Latency = resp_lat) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggbarplot(x = "Age", y = "Response_Latency", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Age",
            palette = c(wes_palette("Zissou1", 10, type = "continuous")[9], "#9EBDBD"),
            #palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Response Latency (sec)" 
       #, subtitle = "STI"
  ) + 
  scale_y_continuous(limits = c(0, 210), expand = c(0, 0))

## NUMBER OF SONGS #####
#plot t test for juv sti vs control - number of songs bar graph
juv_n_song <- 
  data %>% filter(age == "juv") %>% 
  dplyr::select(Condition = condition, Songs = songs) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>% 
  ggbarplot(x = "Condition", y = "Songs",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Number of Songs" 
       #,subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0))  +
  geom_signif(comparisons = list(c("Control", "STI")), annotations = "***",
              y_position = 89, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

#plot t test for STI adults vs juveniles - number of songs bar graph
sti_n_song <- 
  data %>% filter(condition == "sti") %>% 
  dplyr::select(Age = age, Songs = songs) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>%
  ggbarplot(x = "Age", y = "Songs", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Age",
            #palette = c("#505053", "#9EBDBD"),
            palette = c(wes_palette("Zissou1", 10, type = "continuous")[9], "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Number of Songs", 
       #subtitle = ""
  ) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0))

## NUMBER OF FLIGHTS ####
#plot t test for juv sti vs control - number of flights bargraph
juv_fl <- 
  data %>% filter(age == "juv") %>% 
  dplyr::select(Condition = condition, Flights = flights) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>% 
  ggbarplot(x = "Condition", y = "Flights",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Number of Flights", 
       #subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 65), expand = c(0, 0)) +
  geom_signif(comparisons = list(c("Control", "STI")), annotations = "****",
              y_position = 57, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

#plot t test for STI adults vs juveniles - number of flights bargraph
sti_fl <- 
  data %>% filter(condition == "sti") %>% 
  dplyr::select(Age = age, Flights = flights) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>% 
  ggbarplot(x = "Age", y = "Flights", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Age",
            palette = c(wes_palette("Zissou1", 10, type = "continuous")[9], "#9EBDBD"),
            #palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Number of Flights", 
       #subtitle = "STI"
  ) +
  scale_y_continuous(limits = c(0, 65), expand = c(0, 0))

## TIME IN 5 METERS ####
#plot t test for juv sti vs control - time in 5M bar graph
juv_t <- 
  data %>% filter(age=="juv") %>% 
  dplyr::select(Condition = condition, Time_in_5M = time) %>% 
  mutate(Condition = dplyr::recode(Condition, sti = "STI", control = "Control")) %>%  
  ggbarplot(x = "Condition", y = "Time_in_5M",
            order = c("Control", "STI"),
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Condition",
            #palette = c("#bfbfbf", "#505053"),
            palette = c("white", wes_palette("Zissou1", 10, type = "continuous")[9]),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Time Spent in 5m (sec)", 
       #subtitle = "Juvenile"
  ) +
  scale_y_continuous(limits = c(0, 660), expand = c(0, 0)) +
  geom_signif(comparisons = list(c("Control", "STI")), annotations = "****",
              y_position = 600, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

#plot t test for STI adults vs juveniles - time in 5M bar graph
sti_t <- 
  data %>% filter(condition=="sti") %>% 
  dplyr::select(Age = age, Time_in_5M = time) %>% 
  mutate(Age = dplyr::recode(Age, juv = "Juvenile", adult = "Adult")) %>%   
  ggbarplot(x = "Age", y = "Time_in_5M", 
            add = c("mean_se", "dotplot"),
            add.params = list(size = .35, alpha = 1, fill = "black"), 
            fill = "Age",
            #palette = c("#505053", "#9EBDBD"),
            palette = c(wes_palette("Zissou1", 10, type = "continuous")[9], "#9EBDBD"),
            position = position_dodge(1)) +
  theme(legend.position = "none", 
        text = element_text(size = 22)) + 
  labs(x = "", 
       y = "Time Spent in 5m (sec)", 
       #subtitle = "STI"
  ) +
  scale_y_continuous(limits = c(0, 660), expand = c(0, 0)) +
  geom_signif(comparisons = list(c("Juvenile", "Adult")), annotations = "*",
              y_position = 600, tip_length = 10, vjust = .5, 
              size = 1, textsize = 10)

## panel plots ####
juv_resp_lat + juv_t + juv_n_song + juv_fl + 
  plot_layout(design = "AB
  						CD", 
              guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")
ggsave(paste0(path4figs, "behav_juv.pdf"), 
       width = 15, height = 15, dpi = 300)


sti_resp_lat + sti_t + sti_n_song + sti_fl +
  plot_layout(design = "AB
  						CD", 
              guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_sep = " ")
ggsave(paste0(path4figs, "behav_sti.pdf"), 
       width = 15, height = 15, dpi = 300)
