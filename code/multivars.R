library(tidyverse)

# load data ####
beh <- read.csv("rev2/juv_beh_05-05-22.csv", header = T) 
head(beh)

conc <- read.csv("rev2/sg_ms_brain.csv", header = T) 
head(conc)

# pca ####

beh_pca <- beh %>% 
  #filter(condition == "sti") %>% 
  select(subject, call_lat:time) %>% 
  column_to_rownames("subject") %>% mutate_all(log1p) %>% 
  prcomp(scale = T)

summary(beh_pca)
biplot(beh_pca)

factoextra::fviz_contrib(beh_pca, choice = "var", axes = 1, top = 10)

# ensemble datasets
aggr_cons <- 
left_join(
  conc %>% filter(steroid %in% c("p4", "b", "dhc")) %>% 
  pivot_wider(names_from = steroid, values_from = concentration), 
  beh_pca$x %>% data.frame() %>% dplyr::select(PC1) %>% 
  rownames_to_column("subject") %>% 
    mutate(subject = tolower(subject)),
  by = c("subject")
)

aggr_cons %>% filter(is.na(PC1)) # check for unmatched values

## lm plots ####
library(wesanderson)
# pc1 ~ corticosterone
aggr_cons %>% 
  mutate(matrix = dplyr::recode(matrix, 
                                ah = "AH", bnst = "BnST", cb = "Cb", 
                                cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", 
                                plasma = "Plasma", poa = "POA", tna = "TnA", 
                                vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(matrix = fct_relevel(matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>% 
  filter(matrix != "NAc") %>% 
  mutate(across('condition', str_replace, 'con', 'Control'),
         across('condition', str_replace, 'sti', 'STI')) %>% 
  #mutate(across(matrix, 
  #              ~factor(., levels = c("ah","bnst","cb")))) %>% 
  ggplot(aes(x = b, y = PC1, shape = condition, color = condition)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() + 
  labs(x = "Corticosterone", y = "Aggression score (PC1)", 
       shape = "Condition", color = "Condition") +
  scale_color_manual(values = c("gray", wes_palette("Zissou1", 10, type = "continuous")[9]),) +
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  facet_wrap(~matrix, scales = "free")

ggsave(paste0(path4figs, "pc1_b_panel.pdf"), bg = "white",
       width = 10, height = 10, dpi = 300)

# pc1 ~  progesterone
aggr_cons %>% 
  mutate(matrix = dplyr::recode(matrix, 
                                ah = "AH", bnst = "BnST", cb = "Cb", 
                                cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", 
                                plasma = "Plasma", poa = "POA", tna = "TnA", 
                                vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(matrix = fct_relevel(matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>% 
  filter(matrix != "NAc") %>% 
  mutate(across('condition', str_replace, 'con', 'Control'),
         across('condition', str_replace, 'sti', 'STI')) %>% 
  ggplot(aes(x = p4, y = PC1, shape = condition, color = condition)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() + 
  labs(x = "Progesterone", y = "Aggression score (PC1)", 
       shape = "Condition", color = "Condition") +
  scale_color_manual(values = c("gray", wes_palette("Zissou1", 10, type = "continuous")[9]),) +
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  facet_wrap(~matrix, scales = "free")

ggsave(paste0(path4figs, "pc1_p4_panel.pdf"), bg = "white",
       width = 10, height = 10, dpi = 300)

# pc1 - dehydrocorticosterone
aggr_cons %>% 
  mutate(matrix = dplyr::recode(matrix, 
                                ah = "AH", bnst = "BnST", cb = "Cb", 
                                cg = "CG", ls = "LS", nac = "NAc", ncm = "NCM", 
                                plasma = "Plasma", poa = "POA", tna = "TnA", 
                                vmh = "VMH", vta = "VTA", wb = "Blood")) %>% 
  mutate(matrix = fct_relevel(matrix, 
                              "Blood", "Plasma", "NAc", 
                              "POA", "AH", "LS", 
                              "BnST", "VMH", "VTA", "CG", "NCM", "TnA", "Cb")) %>% 
  filter(matrix != "NAc") %>% 
  mutate(across('condition', str_replace, 'con', 'Control'),
         across('condition', str_replace, 'sti', 'STI')) %>% 
  #mutate(across(matrix, 
  #              ~factor(., levels = c("ah","bnst","cb")))) %>% 
  ggplot(aes(x = dhc, y = PC1, shape = condition, color = condition)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() + 
  labs(x = "Dehydrocorticosterone", y = "Aggression score (PC1)", 
       shape = "Condition", color = "Condition") +
  scale_color_manual(values = c("gray", wes_palette("Zissou1", 10, type = "continuous")[9]),) +
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  facet_wrap(~matrix, scales = "free")

ggsave(paste0(path4figs, "pc1_dhc_panel.pdf"), bg = "white",
       width = 10, height = 10, dpi = 300)

## lm regression estimates to csv ####
require(broom)

# corticosterone
aggr_cons %>% 
  filter(condition == "sti") %>% 
  group_by(matrix) %>% 
  do(tidy(lm(data = ., formula = PC1 ~ b*condition))) %>% 
  data.frame() %>% 
  write.csv("pca1_b_lm_pars.csv", row.names = F)

# progesterone
aggr_cons %>% 
  filter(condition == "sti", 
         !matrix %in% c("bnst", "cb", "nac")) %>%
  group_by(matrix) %>% 
  do(tidy(lm(data = ., formula = PC1 ~ p4))) %>% 
  data.frame() %>% mutate(across(where(is.numeric), round, 3)) %>% 
  write.csv("pca1_p4_lm_pars.csv", row.names = F)

# dehydrocorticosterone
aggr_cons %>% 
  filter(condition == "sti",
         !matrix %in% "nac") %>% 
  group_by(matrix) %>% 
  do(tidy(lm(data = ., formula = PC1 ~ dhc))) %>% 
  data.frame() %>% mutate(across(where(is.numeric), round, 3)) %>% 
  write.csv("pca1_dhc_lm_pars.csv", row.names = F)

# wilcox-test on pc1 pc2 with respect to rater ####
left_join(beh_pca$x %>% data.frame() %>% 
            dplyr::select(PC1, PC2) %>% 
            rownames_to_column("subject") %>% 
            mutate(subject = tolower(subject)), 
          beh %>% mutate(subject = tolower(subject)), 
          by = c("subject")) %>% 
  dplyr::select(PC1, PC2, rater) %>% 
  #select_if(is.numeric) %>% sapply(shapiro.test)
  gather(key, value, -rater) %>% 
  group_by(key) %>% 
  do(tidy(wilcox.test(x = .$value, g = .$rater, exact = F))) %>% 
  dplyr::select(-c(method,alternative))