library(tidyverse)
library(reshape2)
library(gplots)
library(RColorBrewer)

########## READ DATA ##########

knownclusterblast <- read_delim("PATH/knownclusterblast.tsv", delim = "\t", col_names = FALSE)
col_names <- c("MAG", "Genome", "Region", "BGC_accession", "BGC_name", "BGC_type")
colnames(knownclusterblast) <- col_names

genome <- read_delim("PATH/genome_information.tsv", delim = "\t") %>% 
  select(MAG, Genome = genome, sample_id = sample, species) %>% 
  mutate(species = ifelse(is.na(species), MAG, species))


########## KNOWN NRPS AND PKS GENE CLUSTERS ##########

# manually edit the names of KnownClusterBlast-identified BGCs
nrps_pks <- knownclusterblast[grepl("(nrp|polyketide)", tolower(knownclusterblast$BGC_type)), ] %>% 
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3}", "\\1", BGC_name)) %>%
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3}", "\\1", BGC_name)) %>%
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3}", "\\1", BGC_name)) %>%
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3}", "\\1", BGC_name)) %>%
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3} / [^ ]+ [A-Z0-9]{1,3}", "\\1", BGC_name)) %>%
  mutate(BGC_name = gsub("^([^ ]+) [A-Z0-9]{1,3}$", "\\1", BGC_name)) %>% 
  mutate(BGC_name = gsub("microcystin-LR / microcystin-FR / microcystin-LA / microcystin-LAba / microcystin-LM / microcystin-LV / microcystin-LL", "microcystin", BGC_name)) %>% 
  distinct(MAG, Genome, BGC_name) # presence/absence of BGCs, the same BGC can only occur once in a genome

num_genome_per_mag <- genome %>% 
  group_by(MAG, species) %>% 
  summarise(number_of_genomes = n()) %>% 
  filter(number_of_genomes > 10) # discard MAGs containing less than 10 genomes

mean_count_gene_clusters <- nrps_pks %>% 
  group_by(MAG, BGC_name) %>% 
  summarise(count = n()) %>% 
  merge(num_genome_per_mag, by = "MAG") %>% 
  mutate(mean_count = count/number_of_genomes)

selected <- mean_count_gene_clusters %>% 
  group_by(BGC_name) %>% 
  summarise(sum = sum(mean_count)) %>% 
  arrange(desc(sum)) %>% 
  head(60) # top 60 most abundant BGC

known_bgc <- merge(mean_count_gene_clusters, selected, by = "BGC_name") %>% 
  mutate(BGC_name= ifelse(BGC_name == "tilivalline / 9-deoxy tilivalline / dihydroxy tilivalline / dehydro tilivalline / dihydroxy-dehydro tilivalline", "dihydroxy tilivalline", BGC_name)) %>% 
  mutate(BGC_name= ifelse(BGC_name == "rosamicin / salinipyrone A / pacificanone A", "rosamicin", BGC_name)) %>% 
  mutate(BGC_name= ifelse(BGC_name == "10,11-dihydro-8-deoxy-12,13-deepoxy-12,13-dihydrochalcomycin", "dihydrochalcomycin", BGC_name))

# convert data frame to matrix
matrix <- dcast(final, species ~ BGC_name, value.var = "mean_count", fill = 0) %>% 
  `row.names<-`(.$species) %>%
  select(-species)
matrix <- data.matrix(matrix)

# heatmap
my_colors <- colorRampPalette(brewer.pal(8, "Reds"))(100)

png("PATH/knownclusterblast.png", width = 400, height = 400, units = "mm", res = 100)
heatmap.2(filtered_mat, dendrogram='none', Rowv=TRUE, Colv=TRUE, trace= 'none', col = my_colors,
          margins = c(10, 15), cexRow = 1, cexCol = 1, key = FALSE)
dev.off()
