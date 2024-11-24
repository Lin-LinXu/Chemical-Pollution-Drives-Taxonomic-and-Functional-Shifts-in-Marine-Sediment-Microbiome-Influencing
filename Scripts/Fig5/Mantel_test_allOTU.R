# Title: Mantel test Fig 5
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

### Loading required package ---------------------------------------------------
library(vegan)
library(ggplot2)
library(reshape)
library(dplyr)
library(vegan)
library(ggpubr)
library(patchwork)
library(rstatix)
library(numbers)
library(tidyverse)
library(ggprism)
library(phyloseq)
### Setting color parameters ---------------------------------------------------
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
#### Microbiome data -----------------------------------------------------------
data_path <- "~/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
Renormalize_data <- function(dataframe) {
  # dataframe <- Procrustes_sign$Profile1[[1]]
  
  otu_tab <- dataframe %>%
    dplyr::rename(Feature = 1) %>%
    column_to_rownames("Feature")
  GPr  = transform_sample_counts(otu_table(otu_tab,taxa_are_rows = TRUE), function(x) x / sum(x) )
  OTU_tab <- data.frame(otu_table(GPr)) %>%
    rownames_to_column("Feature")
  
  return(OTU_tab)
}   
read_taxon <- function(data){
  
  taxon_tmp <- read_tsv(data)
  taxon <- Renormalize_data(taxon_tmp) %>%
    filter(!Feature %in% "Unclassified")
  
  taxon1 <- taxon %>%
    column_to_rownames(var = "Feature") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
  
  taxon2 <- as.matrix(taxon1)
  
  return(taxon2)
}
genus_df <- read_taxon(paste0(data_path,"genus.otu")) 
species_df <- read_taxon(paste0(data_path,"species.otu"))
### Function data --------------------------------------------------------------
Functional_list <- readRDS("~/Desktop/marine_sediment/results/Functional_beta_diversity.rds")
### Small animal  ------------------------------------------------------------
Animal_OTU_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/OTU_2019.csv")
colnames(Animal_OTU_file)[1] <- "ID"
Animal_OTU <- Animal_OTU_file %>%
  column_to_rownames("ID") %>%
  type_convert()
Animal_OTU[Animal_OTU > 0] <- 1
Animal_OTU <- Animal_OTU %>%
  rownames_to_column("ID")
Animal_SAMP <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/SAMP_Linlin.csv")
colnames(Animal_SAMP)[1] <- "ID"
dim(Animal_SAMP)
Animal_SAMP <- Animal_SAMP %>%
  dplyr::select(.,-ARMS_ID,-Ret_ID)
Animal_TAXA_file <- read_tsv("~/Desktop/marine_sediment/data/Shelby_data/TAXA_2019_onlysediment.tsv")
Animal_TAXA <- Animal_TAXA_file %>%
  dplyr::select(.,-ID,-indtype,-SourceID) %>%
  dplyr::rename(ID = Seq) %>%
  relocate(ID,.before = Kingdom) %>%
  mutate(OTU = Species) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Genus),paste0(Genus," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Family),paste0(Family," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Order),paste0(Order," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Class),paste0(Class," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Phylum),paste0(Phylum," OTU",ID),OTU))
#colnames(Animal_TAXA)[1] <- "ID"
dim(Animal_TAXA)
Animal_OTU_df <- Animal_OTU %>%
  left_join(Animal_SAMP) %>%
  relocate(c(Station_ID,Month_Ret), .after = ID) %>%
  mutate(Month_Ret = str_remove_all(Month_Ret,"[INUM]")) %>%
  mutate(ID = str_remove_all(ID,"ARMS_")) %>%
  unite("SampleID",c("Station_ID","Month_Ret","ID"))

Animal_TAXA_tmp <- Animal_TAXA %>%
  dplyr::select(ID,OTU) %>%
  dplyr::rename(Feature = ID)

Animal_otu <- Animal_OTU_df %>%
  filter(SampleID %in% rownames(species_df)) %>%
  #filter(SampleID %in% rownames(COG_path)) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  filter(Feature %in% Animal_TAXA$ID) %>%
  as.data.frame() %>%
  left_join(Animal_TAXA_tmp) %>%
  dplyr::select(.,-Feature) %>%
  dplyr::rename(Feature = OTU) %>%
  relocate(Feature,.before = 1) %>%
  group_by(Feature) %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE)) %>%
  column_to_rownames("Feature") %>%
  t()
## Mantel test  ----------------------------------------------------------------
library(vegan)
spe.dist <- vegdist(species_df, method = 'bray')
genus.dist <- vegdist(genus_df, method = 'bray')
animal.dist <- vegdist(Animal_otu, method = 'jaccard')
cog.dist <- Functional_list$LogEuclidean_dist[[1]]
ko.dist <- Functional_list$LogEuclidean_dist[[2]]

all(colnames(as.matrix(spe.dist)) == colnames(as.matrix(animal.dist)))

spe_mantel<- vegan::mantel(spe.dist, animal.dist, method = 'spearman', permutations = 1000)
genus_mantel<- vegan::mantel(genus.dist, animal.dist, method = 'spearman', permutations = 1000)
cog_mantel<- vegan::mantel(cog.dist, animal.dist, method = 'spearman', permutations = 1000)
ko_mantel<- vegan::mantel(ko.dist, animal.dist, method = 'spearman', permutations = 1000)

spe.label <- paste("Mantel test r:",round(spe_mantel$statistic,digits = 2),
                   ", p = ",signif(spe_mantel$signif,digits = 2))
genus.label <- paste("Mantel test r:",round(genus_mantel$statistic,digits = 2),
                   ", p = ",signif(genus_mantel$signif,digits = 2))
cog.label <- paste("Mantel test r:",round(cog_mantel$statistic,digits = 2),
                     ", p = ",signif(cog_mantel$signif,digits = 2))
ko.label <- paste("Mantel test r:",round(ko_mantel$statistic,digits = 2),
                   ", p = ",signif(ko_mantel$signif,digits = 2))
# ecodist::mantel(animal.dist~spe.dist, mrank =T, nperm = 1000)
# ecodist::mantel(animal.dist~cog.dist, mrank =T, nperm = 1000)
# ecodist::mantel(animal.dist~cog.dist, mrank =T, nperm = 1000)
Plot_df <- data.frame(Species = as.vector(spe.dist),
                      Cog = as.vector(cog.dist),
                      KO = as.vector(ko.dist),
                      Animal = as.vector(animal.dist),
                      Genus = as.vector(genus.dist))
P1 <- ggplot(Plot_df,aes(x = Species, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= spe.label,hjust = 1, vjust = 1,
           size = 16/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("Species Bray-Curtis Dissimilarity") +
  theme_prism(base_size = 16)
P2 <- ggplot(Plot_df,aes(x = Genus, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= genus.label,hjust = 1, vjust = 1,
           size = 16/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("Genus Bray-Curtis Dissimilarity") +
  theme_prism(base_size = 16)
P3 <- ggplot(Plot_df,aes(x = Cog, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= cog.label,hjust = 1, vjust = 1,
           size = 16/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("COG Euclidean Dissimilarity") +
  theme_prism(base_size = 16)
P4 <- ggplot(Plot_df,aes(x = KO, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= cog.label,hjust = 1, vjust = 1,
           size = 16/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("KO Euclidean Dissimilarity") +
  theme_prism(base_size = 16)
p <- P2 + P1 + P3  + 
  plot_layout(ncol = 3) +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 
pdf("~/Desktop/marine_sediment/results/Benthos/Compare_Dissimilarity_all_OTU.pdf",width=14, height= 5)
p
dev.off()

p1 <- P2 + P1 + P3  + P4 +
  plot_layout(ncol = 4) +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 
pdf("~/Desktop/marine_sediment/results/Benthos/Compare_Dissimilarity_withKO.pdf",width=18, height= 5)
p1
dev.off()

