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
data_path <- "~/Desktop/marine_sediment/data/"
### Setting color parameters ---------------------------------------------------
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
#### Microbiome data -----------------------------------------------------------------
read_taxon <- function(data){
  taxon <- read_tsv(data) %>%
    dplyr::select(.,-PC_S_51,-PC_W_33) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as.tibble() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.tibble() %>%
    mutate_at(c(2:length(.)), as.numeric) %>%
    arrange(SampleID) %>%
    column_to_rownames("SampleID")
  return(taxon)
}
genus_df <- read_taxon(paste0(data_path,"genus.otu")) 
species_df <- read_taxon(paste0(data_path,"species.otu"))
### COG data --------------------------------------------------------------
# Assembly prok
COG_file_assembly_prok <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_cog.flt.hellinger.tsv"
COG_path_file_assembly_prok <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_pathway_stat.flt.hellinger.tsv"
COG_assembly_prok <- read_tsv(COG_file_assembly_prok) %>%
  arrange(rowname) %>%
  column_to_rownames()
COG_path_assembly_prok <- read_tsv(COG_path_file_assembly_prok) %>%
  column_to_rownames()
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
#Animal_TAXA <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/TAXA_Linlin.csv")
Animal_TAXA_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/Taxa_2019.csv")
Animal_TAXA <- Animal_TAXA_file %>%
  filter(SourceID %in% "95_Match") %>%
  dplyr::select(.,-`...1`,-indtype,-SourceID) %>%
  rename(ID = Seq) %>%
  relocate(ID,.before = Kingdom)
#colnames(Animal_TAXA)[1] <- "ID"
dim(Animal_TAXA)
Animal_TAXA_df <- Animal_TAXA %>%
  unite(label,Kingdom:Species)
Animal_OTU_df <- Animal_OTU %>%
  left_join(Animal_SAMP) %>%
  relocate(c(Station_ID,Month_Ret), .after = ID) %>%
  mutate(Month_Ret = str_remove_all(Month_Ret,"[INUM]")) %>%
  mutate(ID = str_remove_all(ID,"ARMS_")) %>%
  unite("SampleID",c("Station_ID","Month_Ret","ID"))
Animal_otu <- Animal_OTU_df %>%
  #filter(SampleID %in% rownames(COG_path)) %>%
  column_to_rownames("SampleID") 
Animal_species_otu <- Animal_OTU_df %>%
  filter(SampleID %in% rownames(COG_assembly_prok)) %>%
  arrange(SampleID) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(rowname %in% Animal_TAXA$ID) %>%
  # left_join(Animal_TAXA,by = c("rowname"= "ID")) %>%
  # filter(!is.na(Species)) %>%
  # group_by(Species) %>%
  # summarise_if(is.double,max) %>%
  # column_to_rownames("Species") %>%
  column_to_rownames() %>%
  t()  %>%
  as.data.frame()
Animal_species_otu_tab <- Animal_species_otu %>%
  rownames_to_column("ID")
write_tsv(Animal_species_otu_tab,"~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/OTU.tsv")
## Mantel test  ----------------------------------------------------------------
#library(ecodist)
library(vegan)
spe.dist <- vegdist(species_df, method = 'bray')
#spe.dist1 <- distance(species_df, method = 'bray')
genus.dist <- vegdist(genus_df, method = 'bray')
animal.dist <- vegdist(Animal_species_otu, method = 'jaccard')
#animal.dist1 <- distance(Animal_species_otu, method = 'jaccard')
cog.dist <- vegdist(COG_assembly_prok, method = 'bray')
spe_mantel<- vegan::mantel(spe.dist, animal.dist, method = 'spearman', permutations = 1000)
genus_mantel<- vegan::mantel(genus.dist, animal.dist, method = 'spearman', permutations = 1000)
cog_mantel<- vegan::mantel(cog.dist, animal.dist, method = 'spearman', permutations = 1000)
spe.label <- paste("Mantel test r:",round(spe_mantel$statistic,digits = 2),
                   ", p = ",signif(spe_mantel$signif,digits = 2))
genus.label <- paste("Mantel test r:",round(genus_mantel$statistic,digits = 2),
                   ", p = ",signif(genus_mantel$signif,digits = 2))
cog.label <- paste("Mantel test r:",round(cog_mantel$statistic,digits = 2),
                     ", p = ",signif(cog_mantel$signif,digits = 2))
# ecodist::mantel(animal.dist~spe.dist, mrank =T, nperm = 1000)
# ecodist::mantel(animal.dist~cog.dist, mrank =T, nperm = 1000)
# ecodist::mantel(animal.dist~cog.dist, mrank =T, nperm = 1000)
Plot_df <- data.frame(Species = as.vector(spe.dist),
                      Cog = as.vector(cog.dist),
                      Animal = as.vector(animal.dist),
                      Genus = as.vector(genus.dist))
P1 <- ggplot(Plot_df,aes(x = Species, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= spe.label,hjust = 1, vjust = 1,
           size = 12/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("Species Bray-Curtis Dissimilarity") +
  theme_prism(base_size = 16)
P2 <- ggplot(Plot_df,aes(x = Genus, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= genus.label,hjust = 1, vjust = 1,
           size = 12/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("Genus Bray-Curtis Dissimilarity") +
  theme_prism(base_size = 16)
P3 <- ggplot(Plot_df,aes(x = Cog, y = Animal)) +
  geom_point() +
  geom_rug() +
  geom_smooth(method = lm) +
  #ggpubr::stat_cor(label.x.npc = "left") +
  annotate("text", x = Inf, y = Inf, label= cog.label,hjust = 1, vjust = 1,
           size = 12/.pt) + 
  ylab("Benthos Jaccard Dissimilarity") +
  xlab("COG Bray-Curtis Dissimilarity") +
  theme_prism(base_size = 16)
p <- P1 + P2 + P3  + 
  plot_layout(ncol = 3) +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 
pdf("~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/Compare_Dissimilarity_1.pdf",width=14, height= 4)
p
dev.off()

