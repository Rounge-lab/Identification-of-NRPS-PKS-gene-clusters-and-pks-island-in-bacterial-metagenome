library(tidyverse)
library("RColorBrewer")

########## DATA ########## 

cluster <- read_delim("PATH/antismash_cluster.tsv", delim = "\t") %>%
  mutate(Gene_ID = str_replace_all(as.character(Gene_ID), "['\\[\\]]", "") %>% strsplit(", ")) %>% 
  mutate(Contig = str_extract(Region, "^[^.]+"))

gene <- read_delim("PATH/antismash_gene.tsv", delim = "\t") %>%
  mutate(NRPS_PKS_Modules = ifelse(!is.na(NRPS_PKS_Modules), str_replace_all(as.character(NRPS_PKS_Modules), "['\\[\\]]", ""), "")) %>% 
  mutate(Amino_Acid = ifelse(!is.na(Amino_Acid), str_replace_all(as.character(Amino_Acid), "['\\[\\]]", ""), "")) %>% 
  mutate(Contig = str_extract(Region, "^[^.]+"))

annotation <- read_delim("PATH/dram_annotation.tsv", delim = "\t")

genome <- read_delim("PATH/genome_information.tsv", delim = "\t") %>% 
  select(MAG, Genome = genome, sample_id = sample)

participant <- read_delim("PATH/data_by_sample.txt", delim="\t") %>% 
  select(sample_id, deltaker_id, Prøvetype) %>% 
  na.omit() %>% 
  mutate(sample_id = gsub("_", "-", sample_id))

screening <- read_csv2("PATH/241122_screening_proc.csv") %>% 
  select(deltaker_id, screening_result = final_result_cat3_serr)

contig <- read_delim("PATH/contig_length.tsv", delim = "\t")


########## DISTRIBUTION OF DIFFERENT TYPES OF BGC ########## 

distribution <- cluster %>% 
  group_by(Protocluster_Category) %>% 
  distinct(Genome, Protocluster_Category) %>% 
  summarise(count = n())
distribution$Protocluster_Category <- factor(distribution$Protocluster_Category, levels = rev(c("NRPS", "PKS", "RiPP", "terpene", "other")))

ggplot(distribution, aes(x = Protocluster_Category, y = count, fill = Protocluster_Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "Number of MAGs") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

mags <- paste0("MAG", formatC(1:1941, width = 4, flag = "0"))
types <- c("NRPS", "other", "PKS", "RiPP", "terpene")
all_mag_types <- expand.grid(MAG = mags, Protocluster_Category = types) %>% 
  mutate(count = 0)

variation <- cluster %>%
  group_by(MAG, Protocluster_Category) %>%
  summarise(count = n()) %>% 
  merge(all_mag_types %>% select(MAG, Protocluster_Category), by = c("MAG", "Protocluster_Category"), all = TRUE) %>% 
  mutate(count = ifelse(is.na(count), 0, count)) %>% 
  merge(num_genome_per_mag, by = "MAG") %>% 
  mutate(mean_count = count/number_of_genomes)
variation$Protocluster_Category <- factor(variation$Protocluster_Category, levels = rev(c("NRPS", "PKS", "RiPP", "terpene", "other")))

ggplot(variation, aes(x = Protocluster_Category, y = mean_count, fill = Protocluster_Category)) +
  geom_boxplot(width = 0.6) +
  labs(title = "", x = "", y = "Mean Abundance per Representative MAG") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"))


########## MANUAL ANNOTATION ########## 

dram_clb <- data.frame(
  Uniref_Annotation = c("Colibactin biosynthesis LuxR family transcriptional regulator ClbR", 
                        "Colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB", 
                        "Colibactin polyketide synthase ClbC", 
                        "Colibactin biosynthesis aminomalonyl-acyl carrier protein ClbE", 
                        "Colibactin biosynthesis dehydrogenase ClbF", "ClbG protein", 
                        "Colibactin non-ribosomal peptide synthetase ClbH", 
                        "Colibactin polyketide synthase ClbI", 
                        "Colibactin non-ribosomal peptide synthetase ClbJ", 
                        "Colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK", 
                        "Colibactin biosynthesis amidase ClbL", 
                        "Colibactin non-ribosomal peptide synthetase ClbN", 
                        "Colibactin biosynthesis thioesterase ClbQ", 
                        "Colibactin self-protection protein ClbS"),
  DRAM_Annotation = c("ClbR", "ClbB", "ClbC", "ClbE", "ClbF", "ClbG", "ClbH", "ClbI", "ClbJ", "ClbK", "ClbL", "ClbN", "ClbQ", "ClbS")
)

manual_clb <- data.frame(
  Uniref_Annotation = c("4'-phosphopantetheinyl transferase", 
                        "Transcriptional regulator, LuxR family", 
                        "3-hydroxyacyl-coa dehydrogenase", 
                        "Butyryl-CoA dehydrogenase", 
                        "Amino acid adenylation domain-containing protein", 
                        "Drug/sodium antiporter", "Putative peptide synthetase", 
                        "PKS_KS domain-containing protein", 
                        "Beta-lactamase"),
  NRPS_PKS_Modules = c("ACPS", "", "", "", toString(c("Condensation_LCL", "AMP-binding", "PCP", "Heterocyclization", "AMP-binding", "PCP")), 
                       "", toString(c("Condensation_Starter", "AMP-binding", "PP-binding")), "PKS_KS", ""),
  Manual_Annotation = c("ClbA", "ClbR", "ClbD", "ClbF", "ClbJ", "ClbM", "ClbN", "ClbO", "ClbP")
)

clb_genes = c("ClbA", "ClbR", "ClbB", "ClbC", "ClbD", "ClbE", "ClbF", "ClbG", "ClbH", "ClbI", "ClbJ", "ClbK", "ClbL", "ClbM", "ClbN", "ClbO", "ClbP", "ClbQ", "ClbS")

# case insensitive
insensitive <- function(fun = left_join) {
  new_fun <- fun
  body(new_fun) <- substitute({
    by <- dplyr:::common_by(by, x, y)
    tmp_by_x <- paste0("_", by$x, "_")
    tmp_by_y <- paste0("_", by$y, "_")
    for (i in seq_along(by$x)) {
      x[[tmp_by_x[[i]]]] <- tolower(x[[by$x[[i]]]])
      y[[tmp_by_y[[i]]]] <- tolower(y[[by$y[[i]]]])
      y[[by$y[[i]]]] <- NULL
    }
    res <- fun(x, y, list(x = tmp_by_x, y = tmp_by_y))
    res[tmp_by_x] <- list(NULL)
    res
  })
  new_fun
}

incomplete_manual_annotation <- merge(cluster %>% filter(Protocluster_Category == "NRPS" | Protocluster_Category == "PKS") %>% distinct(MAG, Genome),
                           gene %>% select(Genome, Region, Gene_ID, NRPS_PKS_Modules, Amino_Acid),
                           by = "Genome") %>% 
  merge(annotation %>% select(Genome, Gene_ID, Uniref_Annotation),
        by = c("Gene_ID", "Genome")) %>% 
  insensitive(left_join)(dram_clb, 
                         by = "Uniref_Annotation") %>% 
  # manually annotate specific DRAM annotations containing specific domains as clb genes
  insensitive(left_join)(manual_clb, 
                         by = c("NRPS_PKS_Modules", "Uniref_Annotation")) %>% 
  # DRAM + manual annotation
  mutate(Clb_Annotation = case_when(
    !is.na(DRAM_Annotation) ~ DRAM_Annotation,
    !is.na(Manual_Annotation) ~ Manual_Annotation,
    TRUE ~ NA_character_
  ))

manual_annotation <- incomplete_manual_annotation %>% 
  # remove manually annotated clbA, clbO and clbP if no other clb genes are present
  group_by(Genome) %>%
  mutate(Manual_Annotation = ifelse(Manual_Annotation == "ClbA" & (!any(DRAM_Annotation %in% clb_genes) | !any(Manual_Annotation %in% clb_genes)), NA, Manual_Annotation)) %>%
  mutate(Manual_Annotation = ifelse(Manual_Annotation == "ClbO" & (!any(DRAM_Annotation %in% clb_genes) | !any(Manual_Annotation %in% clb_genes)), NA, Manual_Annotation)) %>%
  mutate(Manual_Annotation = ifelse(Manual_Annotation == "ClbP" & (!any(DRAM_Annotation %in% clb_genes) | !any(Manual_Annotation %in% clb_genes)), NA, Manual_Annotation)) %>%
  ungroup() %>%
  # DRAM + manual annotation
  mutate(Clb_Annotation = case_when(
    !is.na(DRAM_Annotation) ~ DRAM_Annotation,
    !is.na(Manual_Annotation) ~ Manual_Annotation,
    TRUE ~ NA_character_
  ))


########## GENE POSITION IN CONTIGS ########## 

# gene position in contigs for DRAM annotated clb genes
clb_modules <- data.frame(
  Gene = c("ClbA", "ClbR", "ClbB", "ClbC", "ClbD", "ClbE", "ClbF", "ClbG", "ClbH", "ClbI", "ClbJ", "ClbK", "ClbL", "ClbM", "ClbN", "ClbO", "ClbP", "ClbQ", "ClbS"),
  NRPS_PKS_Modules = c("ACPS", "", toString(c("Condensation_DCL", "AMP-binding", "PCP", "PKS_KS", "PKS_AT", "PKS_KR", "PKS_DH2", "PKS_ER", "PCP")), 
                       toString(c("PKS_KS", "PP-binding")), "", "PP-binding", "", "PKS_AT", 
                       toString(c("AMP-binding", "Condensation_LCL", "AMP-binding", "PCP")), toString(c("PKS_KS", "PKS_AT", "PCP")), 
                       toString(c("Condensation_LCL", "AMP-binding", "PCP", "Heterocyclization", "AMP-binding", "PCP")), 
                       toString(c("PKS_KS", "PCP", "Heterocyclization", "AMP-binding", "PCP")), "", "", 
                       toString(c("Condensation_Starter", "AMP-binding", "PP-binding")), "PKS_KS", "", "Thioesterase", "")
)

gene_position <- merge(manual_annotation, 
                       cluster %>% distinct(Genome, Region, Original_Start, Original_End, Contig),
                       by = c("Genome", "Region")) %>% 
  merge(gene %>% select(Genome, Region, Gene_ID, Gene_Location_Start, Gene_Location_End),
        by = c("Genome","Region",  "Gene_ID")) %>% 
  merge(contig,
        by = c("MAG", "Genome", "Contig")) %>% 
  mutate(Original_Gene_Start = Original_Start+Gene_Location_Start) %>% 
  mutate(Original_Gene_End = Original_Start+Gene_Location_End) %>% 
  mutate(contig_start_end = abs(Length - Original_Gene_End) <= 3 | Original_Gene_Start <= 3) %>% 
  filter(grepl("clb", DRAM_Annotation, ignore.case = TRUE)) %>% 
  select(MAG, Genome, Gene_ID, NRPS_PKS_Modules, Amino_Acid, Original_Gene_Start, Original_Gene_End, Length, Clb_Annotation, Uniref_Annotation, DRAM_Annotation, contig_start_end)

incomplete_genes <- merge(gene_position, clb_modules, by.x = "Clb_Annotation", by.y = "Gene") %>%
  mutate(NRPS_PKS_Modules.x = case_when(DRAM_Annotation == "ClbC" & NRPS_PKS_Modules.x == "PKS_KS, PCP" ~ "PKS_KS, PP-binding", TRUE ~ NRPS_PKS_Modules.x)) %>% 
  filter(NRPS_PKS_Modules.x != NRPS_PKS_Modules.y) %>% 
  filter(contig_start_end == FALSE)


########## DRAM vs. MANUAL ANNOTATION ########## 

# frequency of DRAM annotated clb genes 
genes <- c("clbA", "clbR", "clbB", "clbC", "clbD", "clbE", "clbF", "clbG", "clbH", "clbI", "clbJ", "clbK", "clbL", "clbM", "clbN", "clbO", "clbP", "clbQ", "clbS")

count_clb <- lapply(genes, function(gene) {
  sum(grepl(gene, manual_annotation$DRAM_Annotation, ignore.case = TRUE))
})
clb_counts <- as.data.frame(do.call(rbind, Map(cbind, Genes = genes, Count = count_clb)))
clb_counts$Genes <- factor(clb_counts$Genes, levels = rev(clb_counts$Genes))
clb_counts$Count <- as.numeric(clb_counts$Count)

ggplot(clb_counts, aes(Genes, Count)) +
  geom_bar(stat = "identity", fill="#F8766D", width = 0.6, color = "black") +
  labs(title = "DRAM Annotation", x = "", y = "Count", fill = "") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

# frequency of clb genes after the first step of manual annotation
count_incomplete_clb <- lapply(genes, function(gene) {
  sum(grepl(gene, incomplete_manual_annotation$Clb_Annotation, ignore.case = TRUE))
})
incomplete_clb_counts <- as.data.frame(do.call(rbind, Map(cbind, Genes = genes, Count = count_incomplete_clb)))
incomplete_clb_counts$Genes <- factor(incomplete_clb_counts$Genes, levels = rev(incomplete_clb_counts$Genes))
incomplete_clb_counts$Count <- as.numeric(incomplete_clb_counts$Count)

ggplot(incomplete_clb_counts, aes(Genes, Count)) +
  geom_bar(stat = "identity", fill="#00BA38", width = 0.6, color = "black") +
  labs(title = "After Initial Manual Annotation", x = "", y = "Count", fill = "") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

# frequency of clb genes after manual annotation was complete
count_final_clb <- lapply(genes, function(gene) {
  sum(grepl(gene, manual_annotation$Clb_Annotation, ignore.case = TRUE))
})
final_clb_counts <- as.data.frame(do.call(rbind, Map(cbind, Genes = genes, Count = count_final_clb)))
final_clb_counts$Genes <- factor(final_clb_counts$Genes, levels = rev(final_clb_counts$Genes))
final_clb_counts$Count <- as.numeric(final_clb_counts$Count)

ggplot(final_clb_counts, aes(Genes, Count, fill = Genes)) +
  geom_bar(stat = "identity", fill="#619CFF", width = 0.6, color = "black") +
  labs(title = "After Complete Manual Annotation", x = "", y = "Count", fill = "") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

# frequency of clb domains after manual annotation was complete
clb_domains <- manual_annotation %>% 
  filter(grepl("clb", Clb_Annotation, ignore.case = TRUE)) %>%
  mutate(NRPS_PKS_Modules = case_when(
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP" ~ "PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP_2",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP, PKS_KS, PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP" ~ "Condensation_DCL, AMP-binding, PCP_1, PKS_KS, PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP_2",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP, PKS_KS, PKS_AT, PKS_KR, PKS_DH2" ~ "Condensation_DCL, AMP-binding, PCP_1, PKS_KS, PKS_AT, PKS_KR, PKS_DH2",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP, PKS_KS, PKS_AT, PKS_KR" ~ "Condensation_DCL, AMP-binding, PCP_1, PKS_KS, PKS_AT, PKS_KR",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP, PKS_KS, PKS_AT" ~ "Condensation_DCL, AMP-binding, PCP_1, PKS_KS, PKS_AT",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP, PKS_KS" ~ "Condensation_DCL, AMP-binding, PCP_1, PKS_KS",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "Condensation_DCL, AMP-binding, PCP" ~ "Condensation_DCL, AMP-binding, PCP_1",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "AMP-binding, PCP, PKS_KS, PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP" ~ "AMP-binding, PCP_1, PKS_KS, PKS_AT, PKS_KR, PKS_DH2, PKS_ER, PCP_2",
    Clb_Annotation == "ClbB" & NRPS_PKS_Modules == "AMP-binding, PCP" ~ "AMP-binding, PCP_1",
    Clb_Annotation == "ClbC" & NRPS_PKS_Modules == "PKS_KS, PCP" ~ "PKS_KS, PP-binding",
    Clb_Annotation == "ClbH" & NRPS_PKS_Modules == "Condensation_LCL, AMP-binding" ~ "Condensation_LCL, AMP-binding_2",
    Clb_Annotation == "ClbH" & NRPS_PKS_Modules == "AMP-binding, PCP" ~ "AMP-binding_2, PCP",
    Clb_Annotation == "ClbH" & NRPS_PKS_Modules == "AMP-binding, Condensation_LCL, AMP-binding, PCP" ~ "AMP-binding_1, Condensation_LCL, AMP-binding_2, PCP",
    Clb_Annotation == "ClbH" & NRPS_PKS_Modules == "AMP-binding, Condensation_LCL, AMP-binding" ~ "AMP-binding_1, Condensation_LCL, AMP-binding_2",
    Clb_Annotation == "ClbH" & NRPS_PKS_Modules == "AMP-binding, Condensation_LCL" ~ "AMP-binding_1, Condensation_LCL",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "PCP, Heterocyclization, AMP-binding, PCP" ~ "PCP_1, Heterocyclization, AMP-binding_2, PCP_2",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "PCP, Heterocyclization, AMP-binding" ~ "PCP_1, Heterocyclization, AMP-binding_2",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "PCP, Heterocyclization" ~ "PCP_1, Heterocyclization",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "Condensation_LCL, AMP-binding, PCP, Heterocyclization, AMP-binding, PCP" ~ "Condensation_LCL, AMP-binding_1, PCP_1, Heterocyclization, AMP-binding_2, PCP_2",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "Condensation_LCL, AMP-binding, PCP, Heterocyclization, AMP-binding" ~ "Condensation_LCL, AMP-binding_1, PCP_1, Heterocyclization, AMP-binding_2",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "Condensation_LCL, AMP-binding, PCP, Heterocyclization" ~ "Condensation_LCL, AMP-binding_1, PCP_1, Heterocyclization",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "Condensation_LCL, AMP-binding" ~ "Condensation_LCL, AMP-binding_1",
    Clb_Annotation == "ClbJ" & NRPS_PKS_Modules == "AMP-binding, PCP, Heterocyclization, AMP-binding, PCP" ~ "AMP-binding_1, PCP_1, Heterocyclization, AMP-binding_2, PCP_2",
    Clb_Annotation == "ClbK" & NRPS_PKS_Modules == "PKS_KS, PCP, Heterocyclization, AMP-binding, PCP" ~ "PKS_KS, PCP_1, Heterocyclization, AMP-binding, PCP_2",
    Clb_Annotation == "ClbK" & NRPS_PKS_Modules == "PKS_KS, PCP, Heterocyclization" ~ "PKS_KS, PCP_1, Heterocyclization",
    Clb_Annotation == "ClbK" & NRPS_PKS_Modules == "PCP, Heterocyclization, AMP-binding, PCP" ~ "PCP_1, Heterocyclization, AMP-binding, PCP_2",
    Clb_Annotation == "ClbN" & NRPS_PKS_Modules == "Condensation_LCL" ~ "Condensation_Starter",
    TRUE ~ NRPS_PKS_Modules
  )) %>% 
  separate_rows(NRPS_PKS_Modules, sep = ", ") %>% 
  mutate(pair = paste(Clb_Annotation, NRPS_PKS_Modules, sep = ", "))

count_domains <- clb_domains %>% 
  group_by(pair) %>% 
  summarise(certain = n()) %>% 
  mutate(uncertain = 0)

count_domains$uncertain[which(count_domains[["pair"]] == "ClbB, PCP_1")] = count_domains$certain[which(count_domains[["pair"]] == "ClbB, PCP")]
count_domains$uncertain[which(count_domains[["pair"]] == "ClbB, PCP_2")] = count_domains$certain[which(count_domains[["pair"]] == "ClbB, PCP")]
count_domains <- count_domains[-which(count_domains[["pair"]] == "ClbB, PCP"), ]
count_domains$uncertain[which(count_domains[["pair"]] == "ClbH, AMP-binding_1")] = count_domains$certain[which(count_domains[["pair"]] == "ClbH, AMP-binding")]
count_domains$uncertain[which(count_domains[["pair"]] == "ClbH, AMP-binding_2")] = count_domains$certain[which(count_domains[["pair"]] == "ClbH, AMP-binding")]
count_domains <- count_domains[-which(count_domains[["pair"]] == "ClbH, AMP-binding"), ]
count_domains$uncertain[which(count_domains[["pair"]] == "ClbK, PCP_1")] = count_domains$certain[which(count_domains[["pair"]] == "ClbK, PCP")]
count_domains$uncertain[which(count_domains[["pair"]] == "ClbK, PCP_2")] = count_domains$certain[which(count_domains[["pair"]] == "ClbK, PCP")]
count_domains <- count_domains[-which(count_domains[["pair"]] == "ClbK, PCP"), ]

count_domains <- count_domains %>% 
  gather("condition", "count", certain:uncertain)

order <- c("ClbA, ACPS", "ClbR, ", "ClbB, Condensation_DCL", "ClbB, AMP-binding", "ClbB, PCP_1", "ClbB, PKS_KS", "ClbB, PKS_AT", "ClbB, PKS_KR", "ClbB, PKS_DH2", "ClbB, PKS_ER", "ClbB, PCP_2", "ClbC, PKS_KS", "ClbC, PP-binding", "ClbD, ", "ClbE, PP-binding", "ClbF, ", "ClbG, PKS_AT", "ClbH, AMP-binding_1", "ClbH, Condensation_LCL", "ClbH, AMP-binding_2", "ClbH, PCP", "ClbI, PKS_KS", "ClbI, PKS_AT", "ClbI, PCP", "ClbJ, Condensation_LCL", "ClbJ, AMP-binding_1", "ClbJ, PCP_1", "ClbJ, Heterocyclization", "ClbJ, AMP-binding_2", "ClbJ, PCP_2", "ClbK, PKS_KS", "ClbK, PCP_1", "ClbK, Heterocyclization", "ClbK, AMP-binding", "ClbK, PCP_2", "ClbL, ", "ClbM, ", "ClbN, Condensation_Starter", "ClbN, AMP-binding", "ClbN, PP-binding", "ClbO, PKS_KS", "ClbP, ", "ClbQ, Thioesterase", "ClbS, ")
count_domains$pair <- factor(count_domains$pair, levels = order)
count_domains$condition <- factor(count_domains$condition, levels = rev(c("certain", "uncertain")))

ggplot(count_domains, aes(x = pair, y = count, fill = condition)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "Count", fill = "") +
  theme_bw() +
  scale_fill_manual(values = c("certain" = "#619CFF", "uncertain" = "darkgrey"),
                    breaks = c("certain", "uncertain")) +
  scale_x_discrete(labels = rep("Clb?", 44)) +
  theme(legend.position = "none",
        text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))


########## ASSOCIATION BETWEEN PKS AND CRC ########## 

genome_participant_screening <- merge(merge(genome, participant, by = "sample_id"), screening, by = "deltaker_id")

clb_per_genome <- manual_annotation %>%
  mutate(num_clb = grepl("clb", Clb_Annotation, ignore.case = TRUE)) %>%
  group_by(MAG, Genome) %>%
  summarize(clb_genes = sum(num_clb, na.rm = TRUE))

pks_positive_genome <- clb_per_genome %>% 
  filter(clb_genes >= 9) %>% 
  mutate(sample_id = sub("(S-\\d+).*", "\\1", Genome))

baseline_samples <- genome_participant_screening %>% 
  mutate(ecoli_positive = (MAG == "MAG1548")) %>% 
  merge(pks_positive_genome, by = c("MAG", "Genome", "sample_id"), all = TRUE) %>% 
  mutate(pks_positive = !is.na(clb_genes)) %>% 
  filter(Prøvetype == "Baseline") %>% 
  select(sample_id, deltaker_id, screening_result, ecoli_positive, pks_positive)

# association between pks+ genomes (from all baseline samples) and crc
baseline_pks_positive <- baseline_samples %>% 
  filter(pks_positive == TRUE) %>% 
  group_by(screening_result) %>% 
  summarise(pks_positive_ecoli = n())

baseline_screening_result <- genome_participant_screening %>% 
  filter(Prøvetype == "Baseline") %>% 
  group_by(screening_result) %>% 
  distinct(deltaker_id) %>% 
  summarise(total = n())

baseline_pks_screening <- merge(baseline_pks_positive, baseline_screening_result, by = "screening_result") %>% 
  mutate(pks_negative = total - pks_positive_ecoli) %>% 
  mutate(percentage = (pks_positive_ecoli/total)*100) %>% 
  select(screening_result, pks_positive_ecoli, pks_negative, total, percentage)
baseline_pks_screening$screening_result <- factor(baseline_pks_screening$screening_result, levels = c("Advanced lesion", "Non-advanced adenoma", "No adenoma"))

ggplot(baseline_pks_screening, aes(x = screening_result, y = percentage, fill = screening_result)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "Prevalence of pks island (%)") +
  geom_text(aes(label = round(percentage, digits = 1)), position = position_stack(vjust = 0.9), size = 2) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

fisher.test(baseline_pks_screening[1:3, 2:3])
chisq.test(baseline_pks_screening[1:3, 2:3])

# association between pks+ genomes (from E.coli containing baseline samples) and crc
baseline_pks_positive_ecoli <- baseline_samples %>% 
  filter(pks_positive == TRUE & ecoli_positive == TRUE) %>% 
  group_by(screening_result) %>% 
  summarise(pks_positive_ecoli = n())

baseline_ecoli_screening_result <- genome_participant_screening %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(MAG == "MAG1548") %>% 
  group_by(screening_result) %>% 
  summarise(total = n())

baseline_pks_screening <- merge(baseline_pks_positive_ecoli, baseline_ecoli_screening_result, by = "screening_result") %>% 
  mutate(pks_negative = total - pks_positive_ecoli) %>% 
  mutate(percentage = (pks_positive_ecoli/total)*100) %>% 
  select(screening_result, pks_positive_ecoli, pks_negative, total, percentage)
baseline_pks_screening$screening_result <- factor(baseline_pks_screening$screening_result, levels = c("Advanced lesion", "Non-advanced adenoma", "No adenoma"))

ggplot(baseline_pks_screening, aes(x = screening_result, y = percentage, fill = screening_result)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  labs(title = "", x = "", y = "Prevalence of pks island (%)") +
  geom_text(aes(label = round(percentage, digits = 1)), position = position_stack(vjust = 0.9), size = 2) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.position = "none")

fisher.test(baseline_pks_screening[1:3, 2:3])
chisq.test(baseline_pks_screening[1:3, 2:3])
