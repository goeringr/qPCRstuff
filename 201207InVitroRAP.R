RAP <- tibble("FF_1" = c(14.26,14.17,14.47,14.28,14.38,13.85,14.03,15.39,11.86,17.72,18.01,11.80),
              "FF_2" = c(14.25,14.13,14.34,14.22,14.32,13.30,13.67,15.09,11.54,17.70,17.89,11.59),
              "FF_3" = c(14.28,14.16,14.59,14.05,14.11,13.15,13.38,15.15,11.38,17.54,17.78,11.49),
              "GAPDH_1" = c(22.63,22.50,22.87,22.56,21.55,21.53,22.06,22.60,26.48,28.02,30.64,26.56),
              "GAPDH_2" = c(22.71,22.40,22.86,22.55,21.46,21.41,21.88,22.42,26.28,28.11,30.83,26.62),
              "GAPDH_3" = c(22.59,22.38,22.75,22.45,21.16,21.26,21.86,22.38,26.11,27.88,30.37,26.43))

RAP <- t(RAP)

colnames(RAP) <- c("Input_WT", "Input_WTctrl", "Input_mut", "Input_mutctrl", "Snat_WT", "Snat_WTctrl", "Snat_mut", "Snat_mutctrl", "Elute_WT", "Elute_WTctrl", "Elute_mut", "Elute_mutctrl")

RAP <- RAP %>% 
  as_tibble(rownames = "gene") %>% 
  separate(gene, into = c("gene", "rep"), sep = "_") %>% 
  gather(-gene, -rep, key = sample, value = Cq) %>% 
  separate(sample, into = c("sample", "ctrl"), sep = "_")

RAP %>% ggplot(aes(x = ctrl, y = -Cq, fill = gene)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample)

RAP %>% spread(gene,Cq) %>%
  mutate(exp = 2**(-(FF-GAPDH))) %>%
  ggplot(aes(x = ctrl, y = exp, fill = ctrl)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample) + 
  labs(y = "2^(FF-GAPDH)")

RAP %>% spread(gene,Cq) %>%
  filter(sample == "Elute") %>% 
  mutate(exp = 2**(-(FF-GAPDH))) %>%
  ggplot(aes(x = ctrl, y = exp, fill = ctrl)) + 
  geom_boxplot() + 
  theme_cowplot() +
  labs(y = "Relative Firefly Enrichment", x = "")  +
  guides(fill = FALSE) 