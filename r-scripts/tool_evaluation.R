library(tidyverse)
library(ggpubr)

########## BGC DETECTION CAPABILITY ##########

bgc <- read.table("PATH/bgc.tsv", sep = "" , header = F , nrows = 100, na.strings ="", stringsAsFactors= F) %>% 
  setNames(c("genome", "antiSMASH", "SanntiS", "GECCO", "DeepBGC")) 

bgc <- bgc %>% 
  pivot_longer(cols = -genome, names_to = "tool", values_to = "count")
bgc$tool <- factor(bgc$tool, levels = c("genome", "antiSMASH", "SanntiS", "GECCO", "DeepBGC"))  

a <- ggplot(bgc, aes(x = tool, y = count, fill = tool)) +
  geom_boxplot(width = 0.6) +
  labs(title = "", x = "", y = "BGC Count", fill = "") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key.size = unit(1, "lines"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"))


########## NRPS AND PKS DETECTION CAPABILITY ##########

nrps_pks_similar <- read.table("PATH/nrps_pks_1.tsv", sep = "" , header = F , nrows = 100, na.strings ="", stringsAsFactors= F) %>% 
  setNames(c("genome", "antiSMASH", "SanntiS", "GECCO", "DeepBGC"))

nrps_pks_similar <- nrps_pks_similar %>% 
  pivot_longer(cols = -genome, names_to = "tool", values_to = "count")
nrps_pks_similar$tool <- factor(nrps_pks_similar$tool, levels = c("antiSMASH", "SanntiS", "GECCO", "DeepBGC"))  
nrps_pks_similar$genome <- factor(nrps_pks_similar$genome, levels = c("GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7", "GA8", "GA9", "GA10"))  

b <- ggplot(nrps_pks_similar, aes(x = genome, y = count, fill = tool)) +
  geom_bar(position="dodge", stat = "identity", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "NRPS & PKS Count", fill = "") +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.5, "lines"),
        text = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"))


nrps_pks_dissimilar <- read.table("PATH/nrps_pks_2.tsv", sep = "" , header = F , nrows = 100, na.strings ="", stringsAsFactors= F) %>% 
  setNames(c("genome", "antiSMASH", "SanntiS", "GECCO", "DeepBGC"))

nrps_pks_dissimilar <- nrps_pks_dissimilar %>% 
  pivot_longer(cols = -genome, names_to = "tool", values_to = "count")
nrps_pks_dissimilar$tool <- factor(nrps_pks_dissimilar$tool, levels = c("antiSMASH", "SanntiS", "GECCO", "DeepBGC"))  
nrps_pks_dissimilar$genome <- factor(nrps_pks_dissimilar$genome, levels = c("GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7", "GA8", "GA9", "GA10"))  

c <- ggplot(nrps_pks_dissimilar, aes(x = genome, y = count, fill = tool)) +
  geom_bar(position="dodge", stat = "identity", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "NRPS & PKS Count", fill = "") +
  #coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text(size = 8, face = "bold"),
        legend.key.size = unit(0.5, "lines"),
        text = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"))

bc <- ggarrange(b, c,
                nrow = 2, ncol = 1,
                labels = c("B", "C"))

ggsave("PATH/tool_evaluation.png",
       plot = ggarrange(a, bc,
                        nrow = 1, ncol = 2,
                        labels = c("A", ""),
                        widths = c(1.5, 2),
                        common.legend = TRUE,
                        legend = "bottom"),
       width = 250, height = 150, units = "mm", dpi = 300)
