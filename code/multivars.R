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

# ensemble datasets
aggr_cons <- 
left_join(
  conc %>% filter(steroid %in% c("p4", "b")) %>% 
  pivot_wider(names_from = steroid, values_from = concentration), 
  beh_pca$x %>% data.frame() %>% select(PC1) %>% 
  rownames_to_column("subject") %>% 
    mutate(subject = tolower(subject)),
  by = c("subject")
)

aggr_cons %>% filter(is.na(PC1)) # check for unmatched values

# lda by rater ####
# plots ####
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
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  facet_wrap(~matrix, scales = "free")

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
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  facet_wrap(~matrix, scales = "free")

# discriminant

# lm regression estimates to csv ####
require(broom)

# corticosterone
aggr_cons %>% 
  filter(condition == "sti") %>% 
  group_by(matrix) %>% 
  do(tidy(lm(data = ., formula = PC1 ~ b))) %>% 
  data.frame() %>% mutate(across(where(is.numeric), round, 3)) %>% 
  write.csv("pca1_b_lm_pars.csv", row.names = F)

# progesterone
aggr_cons %>% 
  filter(condition == "sti") %>% 
  group_by(matrix) %>% 
  do(tidy(lm(data = ., formula = PC1 ~ p4))) %>% 
  data.frame() %>% mutate(across(where(is.numeric), round, 3)) %>% 
  write.csv("pca1_p4_lm_pars.csv", row.names = F)

