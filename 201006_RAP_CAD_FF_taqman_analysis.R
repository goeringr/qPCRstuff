RAP <- tibble("RE_1" = c(24.72973005,	27.41638802,	30.21187709,	29.78270088,	29.93200576,	30.76525385,	29.14438178,	33.64178028),
"RE_2" = c(24.77654548,	27.35481663,	30.30896017,	29.57487769,	30.00144287,	30.80497949,	28.88589797,	33.58402924),
"RE_3" = c(24.79749591,	27.44371327,	30.29457024,	29.74896962,	29.93626266,	30.89818518,	28.95169812,	34.06560063),
"FF_1" = c(24.53689934,	27.02501995,	29.7237198,	30.04110048,	31.25586348,	30.46988371,	24.48362589,	36.03548358),
"FF_2" = c(24.70372937,	26.96894848,	29.71753737,	29.93585434,	31.06374573,	30.47811725,	24.42645223,	35.01295165),
"FF_3" = c(24.61887799,	27.0509675,	29.50666626,	30.00798738,	31.13181301,	30.53447725,	24.33294819,	35.90207975))

RAP <- t(RAP)

colnames(RAP) <- c("Input_1", "Input_2", "Probe_1", "Probe_2", "Flowthru_1", "Flowthru_2", "Elute_1", "Elute_2")

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
  mutate(exp = 2**(-(FF-RE))) %>%
  ggplot(aes(x = ctrl, y = exp, fill = ctrl)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample) + 
  labs(y = "2^(FF-RE)")

RAP %>% spread(gene,Cq) %>%
  filter(sample == "Elute") %>% 
  mutate(exp = 2**(-(FF-RE))) %>%
  ggplot(aes(x = ctrl, y = exp, fill = ctrl)) + 
  geom_boxplot() + 
  theme_cowplot() +
  labs(y = "Relative Firefly Enrichment", x = "") +
  scale_fill_manual(values = c("#94C47D","#F6BF93")) +
  guides(fill = FALSE) +
  scale_x_discrete(labels=c("+ Probe", "- Probe"))


###################### round 2

RAP2 <- tibble("FF_1" = c(25.08, 39.77, 27.07, 36.85, 27.80, 37.22, 24.79, 36.24),
               "FF_2" = c(24.82, 36.86, 27.06, 40, 27.94, 40, 24.73, 38.62),
              "FF_3" = c(24.67, 36.62, 27.17, 36.41, 27.89, 38.06, 24.88, 37.19),
              "RE_1" = c(25.20, 35.30, 28.25, 35.84, 28.31, 34.51, 29.80, 36.31),
              "RE_2" = c(24.85, 34.81, 28.20, 37.19, 28.44, 36.40, 29.68, 37.40),
              "RE_3" = c(24.81, 34.79, 28.42, 35.92, 28.48, 35.61, 29.93, 35.07))

RAP2 <- t(RAP2)

colnames(RAP2) <- c("Input_FF", "Input_LOX", "Probe_FF", "Probe_LOX", "Flowthru_FF", "Flowthru_LOX", "Elute_FF", "Elute_LOX")

RAP2 <- RAP2 %>% 
  as_tibble(rownames = "gene") %>% 
  separate(gene, into = c("gene", "rep"), sep = "_") %>% 
  gather(-gene, -rep, key = sample, value = Cq) %>% 
  separate(sample, into = c("sample", "cell"), sep = "_")

RAP2 %>% ggplot(aes(x = cell, y = -Cq, fill = gene)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample)

RAP2 %>% spread(gene, Cq) %>%
  mutate(exp = 2**(-(FF-RE))) %>%
  ggplot(aes(x = cell, y = exp, fill = cell)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample) + 
  labs(y = "2^(FF-RE)")

RAP2 %>% spread(gene,Cq) %>%
  filter(sample == "Elute") %>% 
  mutate(exp = 2**(-(FF-RE))) %>%
  ggplot(aes(x = cell, y = exp, fill = cell)) + 
  geom_boxplot() + 
  theme_cowplot() +
  labs(y = "Relative Firefly Enrichment", x = "") +
  scale_fill_manual(values = c("#94C47D","#F6BF93")) +
  guides(fill = FALSE) +
  scale_x_discrete(labels=c("Firefly", "no Firefly"))

###################### rep1??

RAP_rep1 <- tibble("RE_1" = c(25.87,25.66,26.87,26.42),
               "RE_2" = c(26.03,25.72,27.13,26.34),
               "RE_3" = c(25.44,NA,26.84,26.33),
               "FF_1" = c(28.62,28.31,26.21,25.44),
               "FF_2" = c(28.77,28.15,26.38,25.53),
               "FF_3" = c(28.91,NA,26.11,25.36))

RAP_rep1 <- t(RAP_rep1) %>% na.omit()

colnames(RAP_rep1) <- c("Input_WT", "Input_Mutant", "Elute_WT", "Elute_Mutant")

RAP_rep1 <- RAP_rep1 %>% 
  as_tibble(rownames = "gene") %>% 
  separate(gene, into = c("gene", "rep"), sep = "_") %>% 
  gather(-gene, -rep, key = sample, value = Cq) %>% 
  separate(sample, into = c("sample", "cell"), sep = "_")

RAP_rep1 %>% ggplot(aes(x = cell, y = -Cq, fill = gene)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample)

RAP_rep1 %>% spread(gene, Cq) %>%
  mutate(exp = 2**(-(FF-RE))) %>%
  ggplot(aes(x = cell, y = exp, fill = cell)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  facet_grid(.~sample) + 
  labs(y = "2^(FF-RE)") +
  scale_fill_manual(values = c("#F6BF93", "#94C47D")) +
  labs(y = "Relative Firefly Enrichment", x = "") +
  scale_x_discrete(limits = c("WT", "Mutant"))

