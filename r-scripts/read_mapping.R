library(tidyverse)
library(ggpubr)

########## DATA ##########

genes <- read_delim("PATH/bedcov_genes.tsv", delim = "\t", col_names = FALSE) %>% 
  setNames(c("clbA", "clbR", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG", "clbH", "clbI", "clbJ", "clbK", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ", "clbS", "sample_id"))

domains <- read_delim("PATH/bedcov_domains.tsv", delim = "\t", col_names = FALSE) %>% 
  setNames(c("clbA,APC", "clbR", "clbB,Condensation_DCL", "clbB,AMP-binding", "clbB,PCP_1", "clbB,PKS_KS", "clbB,PKS_AT", "clbB,PKS_KR", "clbB,PKS_DH2", "clbB,PKS_ER", "clbB,PCP_2", "clbC,PKS_KS", "clbC,PP-binding", "clbD", "clbE,PP-binding", "clbF", "clbG,PKS_AT", "clbH,AMP-binding_1", "clbH,Condensation_LCL", "clbH,AMP-binding_2", "clbH,PCP", "clbI,PKS_KS", "clbI,PKS_AT", "clbI,PCP", "clbJ,Condensation_LCL", "clbJ,AMP-binding_1", "clbJ,PCP_1", "clbJ,Heterocyclization", "clbJ,AMP-binding_2", "clbJ,PCP_2", "clbK,PKS_KS", "clbK,PCP_1", "clbK,Heterocyclization", "clbK,AMP-binding", "clbK,PCP_2", "clbL", "clbM", "clbN,Condensation_Starter", "clbN,AMP-binding_2", "clbN,PP-binding", "clbO,PKS_KS", "clbP", "clbQ,Thioesterase", "clbS", "sample_id"))

gene_length <- read.table("PATH/clb_genes.bed", sep = "" , header = F , nrows = 100, na.strings ="", stringsAsFactors= F) %>% 
  setNames(c("ref_accession", "start_position", "end_position", "gene")) %>% 
  mutate(length = end_position - start_position) %>% 
  select(gene, length)

domain_length <- read.table("PATH/clb_domains.bed", sep = "" , header = F , nrows = 100, na.strings ="", stringsAsFactors= F) %>% 
  setNames(c("ref_accession", "start_position", "end_position", "module")) %>% 
  mutate(length = end_position - start_position) %>% 
  select(module, length)

coverage_pks <- read_delim("PATH/coverage_pks.tsv", delim = "\t", col_names = FALSE) %>% 
  setNames(c("covered_bases", "coverage", "mean_depth_of_coverage", "sample_id"))

coverage_ecoli <- read_delim("PATH/coverage_ecoli.tsv", delim = "\t", col_names = FALSE) %>% 
  setNames(c("covered_bases", "coverage", "mean_depth_of_coverage", "sample_id"))

idxstats_pks <- read_delim("PATH/idxstats_pks.tsv", delim = "\t", col_names = FALSE) %>% 
  setNames(c("ref_seq_length", "mapped_reads", "placed_but_unmapped_reads", "unmanpped_reads", "total_reads", "sample_id")) %>% 
  select(sample_id, mapped_reads)

participant <- read_delim("PATH/data_by_sample.txt", delim="\t") %>% 
  select(sample_id, deltaker_id, Prøvetype, total_reads = Total_Reads_QC_ATLAS) %>% 
  na.omit() %>% 
  mutate(sample_id = gsub("_", "-", sample_id))

screening <- read_csv2("PATH/241122_screening_proc.csv") %>% 
  select(deltaker_id, screening_result = final_result_cat3_serr)


########## READ DEPTH OF CLB GENES ##########

# get the mean depth for each clb gene
gene_mean_depth <- map2_dfc(genes[,1:19], gene_length$length, `/`) %>% 
  cbind(genes[, -(1:19)])

# get the read depth per million for each clb gene
gene_read_depth_per_million <- merge(gene_mean_depth, participant %>% select(sample_id, total_reads), by = "sample_id") %>% 
  mutate_at(vars(2:20), ~./total_reads * 1e6) %>% 
  select(-total_reads) %>% 
  filter(rowSums(select(., 2:20)) != 0)

# convert the data frame format
sample_gene_depth <- gene_read_depth_per_million %>% 
  pivot_longer(cols = -sample_id, names_to = "gene", values_to = "depth")

# order the genes
clb_genes <- c("clbA", "clbR", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG", "clbH", "clbI", "clbJ", "clbK", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ", "clbS")
sample_gene_depth$gene <- factor(sample_gene_depth$gene, levels = clb_genes)

# pks+ samples
pks_positive_samples <- coverage_pks %>% 
  filter(coverage > 50)

sample_gene_depth <- sample_gene_depth %>% 
  inner_join(pks_positive_samples, by = "sample_id") %>%
  select(names(sample_gene_depth)) %>% 
  mutate(depth = depth+0.01) # added with 0.01 to not discard log(0)

# boxplot
ggplot(sample_gene_depth, aes(x = gene, y = depth, fill = gene)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  labs(title = "", x = "", y = "Read Depth per Million") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  scale_y_log10()


########## READ DEPTH OF CLB DOMAINS ##########

# get the mean depth for each clb domain
domain_mean_depth <- map2_dfc(domains[,1:44], domain_length$length, `/`) %>% 
  cbind(domains[, -(1:44)])

# get the read depth per million for each clb domain
domain_read_depth_per_million <- merge(domain_mean_depth, participant %>% select(sample_id, total_reads), by = "sample_id") %>% 
  mutate_at(vars(2:45), ~./total_reads * 1e6) %>% 
  select(-total_reads) %>% 
  filter(rowSums(select(., 2:20)) != 0)

# convert the data frame format
sample_domain_depth <- domain_read_depth_per_million %>% 
  pivot_longer(cols = -sample_id, names_to = "domain", values_to = "depth")

# order the domains
clb_domains <- c("clbA,APC", "clbR", "clbB,Condensation_DCL", "clbB,AMP-binding", "clbB,PCP_1", "clbB,PKS_KS", "clbB,PKS_AT", "clbB,PKS_KR", "clbB,PKS_DH2", "clbB,PKS_ER", "clbB,PCP_2", "clbC,PKS_KS", "clbC,PP-binding", "clbD", "clbE,PP-binding", "clbF", "clbG,PKS_AT", "clbH,AMP-binding_1", "clbH,Condensation_LCL", "clbH,AMP-binding_2", "clbH,PCP", "clbI,PKS_KS", "clbI,PKS_AT", "clbI,PCP", "clbJ,Condensation_LCL", "clbJ,AMP-binding_1", "clbJ,PCP_1", "clbJ,Heterocyclization", "clbJ,AMP-binding_2", "clbJ,PCP_2", "clbK,PKS_KS", "clbK,PCP_1", "clbK,Heterocyclization", "clbK,AMP-binding", "clbK,PCP_2", "clbL", "clbM", "clbN,Condensation_Starter", "clbN,AMP-binding_2", "clbN,PP-binding", "clbO,PKS_KS", "clbP", "clbQ,Thioesterase", "clbS")
sample_domain_depth$domain <- factor(sample_domain_depth$domain, levels = clb_domains)

# pks+ samples
pks_positive_samples <- coverage_pks %>% 
  filter(coverage > 50)

sample_domain_depth <- sample_domain_depth %>% 
  inner_join(pks_positive_samples, by = "sample_id") %>%
  select(names(sample_domain_depth)) %>% 
  mutate(depth = depth+0.01) # added with 0.01 to not discard log(0)

# boxplot
ggplot(sample_domain_depth, aes(x = domain, y = depth, fill = domain)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  labs(title = "", x = "", y = "Read Depth per Million") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  scale_y_log10()


########## IDENTIFICATION OF PKS ISLAND AND E. COLI ##########

# dotplot
pks_ecoli_coverage <- merge(coverage_pks, coverage_ecoli, by = "sample_id")

ggplot(pks_ecoli_coverage, aes(coverage.x, coverage.y)) +
  geom_point(color = "lightblue") +
  labs(title = "", x = "pks island Coverage", y = "E. coli Coverage") +
  theme_bw() +
  scale_color_manual(values = c("highlight" = "salmon")) + 
  geom_point(data = pks_ecoli_coverage %>% filter(coverage.x>50 & coverage.y<50), aes(color = "highlight")) +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black")

# 184 pks+, 657 e.coli, 175 pks+ e.coli
pks_positive <- coverage_pks %>% 
  filter(coverage > 50)
ecoli_positive <- coverage_ecoli %>% 
  filter(coverage > 50)
pks_positive_ecoli <- merge(pks_positive, ecoli_positive, by = "sample_id") # 175


########## ASSOCIATION BETWEEN PKS AND CRC ##########

# correlation between pks+ (from baseline samples) and screening result
pks_crc <- merge(participant %>% select(-total_reads),
                 screening,
                 by = "deltaker_id") %>% 
  merge(pks_positive,
        by = "sample_id") %>% 
  filter(Prøvetype == "Baseline")

all <- merge(merge(participant %>% select(-total_reads), screening, by = "deltaker_id"), genes, by = "sample_id")
screening_result <- all %>% # 1034 baseline samples
  filter(Prøvetype == "Baseline") %>% 
  group_by(screening_result) %>% 
  summarise(total = n())

pks_screening <- merge(pks_crc %>% group_by(screening_result) %>% summarise(pks_positive = n()), 
                       screening_result, 
                       by = "screening_result", all.y = TRUE) %>% 
  mutate(percentage = (pks_positive/total)*100)
pks_screening$screening_result <- factor(pks_screening$screening_result, levels = c("Advanced lesion", "Non-advanced adenoma", "No adenoma"))

ggplot(pks_screening, aes(x = screening_result, y = percentage, fill = screening_result)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "", x = "", y = "prevalence of pks island (%)") +
  geom_text(aes(label = round(percentage, digits = 1)), position = position_stack(vjust = 0.9), size = 2) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

fisher.test(pks_screening[1:3, 2:3])
chisq.test(pks_screening[1:3, 2:3])

# distribution of read count
read_count <- merge(pks_crc %>% select(sample_id, screening_result), 
                    idxstats_pks, by = "sample_id") %>% 
  merge(participant %>%  select(sample_id, total_reads), by = "sample_id") %>% 
  mutate(normalized_reads = (mapped_reads/total_reads)*1e6) %>% 
  select(-mapped_reads, -total_reads)
read_count$screening_result <- factor(read_count$screening_result, levels = c("Advanced lesion", "Non-advanced adenoma", "No adenoma"))

ggplot(read_count, aes(x = screening_result, y = normalized_reads, fill = screening_result)) +
  geom_boxplot(width = 0.6) +
  labs(title = "", x = "", y = "read count") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

kruskal.test(normalized_reads ~ screening_result, data = read_count)
