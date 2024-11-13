## Loading required package ---------------------------------------------------
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
## Setting color parameters ---------------------------------------------------
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                   Gradient = c("High","High","High","High","Low","Low","Low"))
## Load data ------------------------------------------------------------------
microbiome_spe <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/alpha_diversity_species_STAMP.rds")
microbiome_genus <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/alpha_diversity_genus_STAMP.rds")
function_lst <- readRDS("~/Desktop/marine_sediment/results/Function/functional_alpha_diversity.rds")
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
  #filter(SampleID %in% rownames(COG_path)) %>%
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
  summarise(across(where(is.numeric), max, na.rm = TRUE))

Animal_species_richness <- Animal_otu %>%
  column_to_rownames("Feature") %>% 
  t() %>%
  as.data.frame() %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID,total) %>%
  filter(SampleID %in% microbiome_spe$alpha_diversity$SampleID) %>%
  dplyr::rename(Richness = total)

## Plot ----------------------------------------------------------------------
Microbiome_Spe_richness <- microbiome_spe$alpha_diversity %>%
  dplyr::select(SampleID,Richness) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

Microbiome_Genus_richness <- microbiome_genus$alpha_diversity %>%
  dplyr::select(SampleID,Richness) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

Microbiome_COG_richness <- function_lst$COG_alphaD_list$alpha_diversity %>%
  dplyr::select(SampleID,Richness) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

Microbiome_KO_richness <- function_lst$KO_alphaD_list$alpha_diversity %>%
  dplyr::select(SampleID,Richness) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

Benthos_richness <- Animal_species_richness %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID")

all(rownames(Microbiome_COG_richness) == rownames(Benthos_richness))

Plot_df <- data.frame(Microbiome_Spe = Microbiome_Spe_richness$Richness,
                      Microbiome_Genus = Microbiome_Genus_richness$Richness,
                      KO = Microbiome_KO_richness$Richness,
                      COG = Microbiome_COG_richness$Richness,
                      Animal = Benthos_richness$Richness,
                      tmp = rownames(Benthos_richness)) %>%
  `rownames<-`(rownames(Benthos_richness)) %>%
  separate(tmp,c("Site","Season","tmp")) %>%
  dplyr::select(.,-tmp) %>%
  mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
  mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC")))
P1 <- ggplot(Plot_df,aes(x = Microbiome_Spe, y = Animal)) +
  geom_point(aes(color = Site,shape = Season)) +
  geom_rug() +
  geom_smooth(method = lm) +
  ggpubr::stat_cor(size = 16/.pt) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_manual(values = cbPalette) +
  ylab("Benthos Richness") +
  xlab("Microbiome species Richness") +
  theme_prism(base_size = 16)
P2 <- ggplot(Plot_df,aes(x = Microbiome_Genus, y = Animal)) +
  geom_point(aes(color = Site,shape = Season)) +
  geom_rug() +
  geom_smooth(method = lm) +
  ggpubr::stat_cor(size = 16/.pt) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_manual(values = cbPalette) +
  ylab("Benthos Richness") +
  xlab("Microbiome genus Richness") +
  theme_prism(base_size = 16)
P3 <- ggplot(Plot_df,aes(x = COG, y = Animal)) +
  geom_point(aes(color = Site,shape = Season)) +
  geom_rug() +
  geom_smooth(method = lm) +
  ggpubr::stat_cor(size = 16/.pt) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_manual(values = cbPalette) +
  ylab("Benthos Richness") +
  xlab("COG Richness") +
  theme_prism(base_size = 16)

p <- P2 + P1 + P3  + 
  plot_layout(ncol = 3) +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 


pdf("~/Desktop/marine_sediment/results/Benthos/Compare_richness_allOTU.pdf",width=12, height= 5)
p
dev.off()

p2 <- P1 + P3  + 
  plot_layout(ncol = 2) +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 

pdf("~/Desktop/marine_sediment/results/Benthos/Compare_richness_allOTU1.pdf",width=12, height= 5)
p2
dev.off()

