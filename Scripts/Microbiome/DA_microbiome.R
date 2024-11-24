# Title: Differential Abundance microbiome
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024

## Loading required packages ---------------------------------------------------
library(Maaslin2)
library(tidyverse)
library(phyloseq)
#library(MatrixExtra)
## Loading data ----------------------------------------------------------------
data_path <- "~/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
genus_df <- read_tsv(paste0(data_path,"genus.otu"))
species_df <- read_tsv(paste0(data_path,"species.otu"))
family_df <- read_tsv(paste0(data_path,"family.otu"))
phylum_df <- read_tsv(paste0(data_path,"phylum.otu"))
## Prepare df ------------------------------------------------------------------
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
Pre_filter <- function(data)
{
  #data <- species_df
  taxon <- Renormalize_data(data) %>% 
    filter(!Feature %in% "Unclassified") %>% 
    dplyr::rename(SampleID = 1) %>%
    dplyr::filter(!SampleID %in% "uncultured archaeon") %>%
    column_to_rownames("SampleID")
  
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  
  meta_data <- data.frame(name = colnames(taxon)) %>%
    mutate(tmp = name) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    left_join(gradient)
  meta_tmp <- meta_data %>% dplyr::select(name,Site)
  
  flt_phyloseq <- function(df) {
    df <- df %>% column_to_rownames("name") 
    GPfr_stage1 = filter_taxa(otu_table(df,taxa_are_rows = FALSE),
                              function(x) sum(x > 0.001) >= 3, TRUE)
    otu_flt_stage1 <- GPfr_stage1@.Data
  }
  
  OTU_tab <- taxon %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    left_join(meta_tmp) %>%
    group_nest(Site) %>%
    mutate(data = map(data, function(df) flt_phyloseq(df))) %>%
    mutate(species = map(data, function(df) colnames(df)))
  
  Passed_species <- unlist(OTU_tab$species) %>% unique()
  
  final_df <- taxon %>%
    rownames_to_column("Feature") %>%
    filter(Feature %in% Passed_species) %>%
    column_to_rownames("Feature") %>%
    t()
  meta_data <- meta_data %>%
    column_to_rownames("name")
  
  Feature_ID <- data.frame(
    name = colnames(final_df),
    ID = paste0("Taxa", 1:length(colnames(final_df)))
  )
  
  colnames(final_df) <- Feature_ID$ID
  
  return_list <- list(
    taxon = final_df,
    meta_data = meta_data,
    Feature_ID = Feature_ID
  )
}
species_list <- Pre_filter(species_df)
all(rownames(species_list$taxon) == rownames(species_list$meta_data))
genus_list <- Pre_filter(genus_df)
all(rownames(genus_list$taxon) == rownames(genus_list$meta_data))
family_list <- Pre_filter(family_df)
all(rownames(family_list$taxon) == rownames(family_list$meta_data))
phylum_list <- Pre_filter(phylum_df)
all(rownames(phylum_list$taxon) == rownames(phylum_list$meta_data))
## Maaslin2 --------------------------------------------------------------------
func_Maaslin2 <- function(data_list){
  fit_data = Maaslin2(input_data     = data_list$taxon, 
                      input_metadata = data_list$meta_data, 
                      min_prevalence = 0,
                      min_abundance = 0,
                      analysis_method = "CPLM",
                      normalization  = "NONE",
                      transform = "NONE",
                      output         = "demo_output", 
                      fixed_effects  = c("Gradient", "Season"),
                      reference      = c("Gradient,Low"))
  
  directoryPath <- "~/demo_output"
  unlink(directoryPath, recursive = TRUE)
  
  sign_df <- fit_data$results %>%
    filter(metadata %in% "Gradient") %>%
    filter(pval <= 0.05, qval <= 0.20) %>%
    mutate()
  
  re_list <- list(
    results = fit_data$results,
    sign_df = sign_df
  )
  return(re_list)
}
## Results Pollution -----------------------------------------------------------
species_relist <- func_Maaslin2(species_list)
genus_relist <- func_Maaslin2(genus_list)
family_relist <- func_Maaslin2(family_list)
phylum_relist <- func_Maaslin2(phylum_list)

genus_df <- genus_relist$sign_df %>%
  left_join(genus_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

species_df <- species_relist$sign_df %>%
  left_join(species_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

family_df <- family_relist$sign_df %>%
  left_join(family_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

phylum_df <- phylum_relist$sign_df %>%
  left_join(phylum_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

Microbiome_DA_List <- list(
  species_relist = species_relist,
  genus_relist = genus_relist,
  species_list = species_list,
  genus_list = genus_list,
  family_relist = family_relist,
  phylum_relist = phylum_relist,
  family_list = family_list,
  phylum_list = phylum_list
)

saveRDS(Microbiome_DA_List,"~/Desktop/marine_sediment/results/MEGAN6/Microbiome_DA_pollution_List.rds")

## Check Nitrospina ---------------
Microbiome_DA_List <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/Microbiome_DA_pollution_List.rds")

phylum_df <- Microbiome_DA_List$phylum_relist$sign_df %>%
  left_join(Microbiome_DA_List$phylum_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name)) %>%
  mutate(rank = "phylum")

genus_df <- Microbiome_DA_List$genus_relist$sign_df %>%
  left_join(Microbiome_DA_List$genus_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name)) %>%
  mutate(rank = "genus")

species_df <- Microbiome_DA_List$species_relist$sign_df %>%
  left_join(Microbiome_DA_List$species_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name)) %>%
  mutate(rank = "species")

re_df <- bind_rows(list(phylum_df,genus_df,species_df))
library(xlsx)
write.xlsx(re_df, "~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiom/Supplementary_tables_new/Supplementary_Tabs.xlsx",
           sheetName="Table S4",row.names = T,append=TRUE)
