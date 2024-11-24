# Title: Differential Abundance benthos
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Loading required packages ---------------------------------------------------
library(Maaslin2)
library(tidyverse)
library(phyloseq)
#library(MatrixExtra)
## Loading data ----------------------------------------------------------------
Animal_OTU_file <- read_csv("/data/Desktop/marine_sediment/data/Shelby_data/OTU_2019.csv")
colnames(Animal_OTU_file)[1] <- "ID"
Animal_OTU <- Animal_OTU_file %>%
  column_to_rownames("ID") %>%
  type_convert()
Animal_OTU[Animal_OTU > 0] <- 1
Animal_OTU <- Animal_OTU %>%
  rownames_to_column("ID")
Animal_SAMP <- read_csv("/data/Desktop/marine_sediment/data/Shelby_data/SAMP_Linlin.csv")
colnames(Animal_SAMP)[1] <- "ID"
dim(Animal_SAMP)
Animal_SAMP <- Animal_SAMP %>%
  dplyr::select(.,-ARMS_ID,-Ret_ID)
Animal_TAXA_file <- read_tsv("/data/Desktop/marine_sediment/data/Shelby_data/TAXA_2019_onlysediment.tsv")
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
## Prepare df ------------------------------------------------------------------
Pre_filter <- function(data)
{
  #data <- Animal_otu
  
  taxon <- data %>% 
    column_to_rownames("Feature")
  
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low")
                     )
  
  meta_data <- data.frame(name = colnames(taxon)) %>%
    mutate(tmp = name) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    left_join(gradient)
  meta_tmp <- meta_data %>% dplyr::select(name,Gradient)
  
  flt_phyloseq <- function(df) {
    df <- df %>% column_to_rownames("name") 
    GPfr_stage1 = filter_taxa(otu_table(df,taxa_are_rows = FALSE),
                              function(x) sum(x > 0) >= 3, TRUE)
    otu_flt_stage1 <- GPfr_stage1@.Data
  }
  
  OTU_tab <- taxon %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("name") %>%
    left_join(meta_tmp) %>%
    group_nest(Gradient) %>%
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
    ID = paste0("Feature", 1:length(colnames(final_df)))
  )
  
  colnames(final_df) <- Feature_ID$ID
  
  return_list <- list(
    taxon = final_df,
    meta_data = meta_data,
    Feature_ID = Feature_ID
  )
}

Animal_otu_list <- Pre_filter(Animal_otu)
all(rownames(Animal_otu_list$taxon) == rownames(Animal_otu_list$meta_data))
## Maaslin2 --------------------------------------------------------------------
func_Maaslin2 <- function(data_list){
  fit_data = Maaslin2(input_data     = data_list$taxon, 
                      input_metadata = data_list$meta_data, 
                      min_prevalence = 0,
                      min_abundance = 0,
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "NONE",
                      output         = "demo_output", 
                      fixed_effects  = c("Gradient", "Season"),
                      reference      = c("Gradient,Low"))
  
  directoryPath <- "/data/demo_output"
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
Animal_otu_relist <- func_Maaslin2(Animal_otu_list)

Animal_otu_df <- Animal_otu_relist$sign_df %>%
  left_join(Animal_otu_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

Animal_List <- list(
  Animal_otu_relist = Animal_otu_relist,
  Animal_otu_list = Animal_otu_list
)

saveRDS(Animal_List,"/data/Desktop/marine_sediment/results/Benthos/Animal_List.rds")

