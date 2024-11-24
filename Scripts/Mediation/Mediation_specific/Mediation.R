# Title: Mediation individual level
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Final results runned on server.
## Load library ----------------------------------------------------------------
library(tidyverse)
library(mediation)
set.seed(123)
## Load data -------------------------------------------------------------------
Input_data <- readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/Input_data.rds")
Function_data <-  readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/functional_profile.rds")
## Significant features --------------------------------------------------------
microbiome_sign <- readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/Microbiome_DA_pollution_List.rds")
benthos_sign <- readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/Animal_List.rds")
Pollution_path_list <- readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/Function_DA_PATHWAY_pollution_List.rds")
Pollution_term_list <- readRDS("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/data/Function_DA_pollution_List.rds")
microbiome_sign_df <- microbiome_sign$species_relist$sign_df %>%
  dplyr::rename(ID = feature, tmp = name) %>%
  left_join(microbiome_sign$species_list$Feature_ID)
benthos_sign_df <- benthos_sign$Animal_otu_relist$sign_df %>%
  dplyr::rename(ID = feature, tmp = name) %>%
  left_join(benthos_sign$Animal_otu_list$Feature_ID)
## Mediation Analyses  pollution -> microbiome -> animal ----------------------------------------------------------
microbiome_mat <- Input_data$species_df %>%
  as.data.frame() %>%
  #dplyr::select(where(~is.numeric(.) && var(., na.rm = TRUE) != 0)) %>%
  dplyr::select(microbiome_sign_df$name) %>%
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID")
Benthos_mat <- Input_data$Animal_otu %>%
  as.data.frame() %>%
  #dplyr::select(where(~is.numeric(.) && var(., na.rm = TRUE) != 0)) %>%
  dplyr::select(benthos_sign_df$name) %>%
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID")
Pollution_mat <- Input_data$Pollution_df %>%
  as.data.frame() %>%
  dplyr::select(where(~is.numeric(.) && var(., na.rm = TRUE) != 0)) %>%
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID")
all(rownames(Benthos_mat) == rownames(Pollution_mat))


ALl_mat <- cbind(microbiome_mat,Benthos_mat,Pollution_mat)

library(foreach)
library(doParallel)
source("~/Desktop/marine_sediment/scripts/Mediation_specific/Mediation_specific/mediation_functions_new.R")
## pollution -> microbiome -> animal -------------------------------------------
registerDoParallel(22)

re_df <- foreach (i = 1:dim(Pollution_mat)[2],.combine=rbind) %dopar% {
  #i <- 3
  set.seed(i)
  print(i)
  a = colnames(Pollution_mat)[i]
  mediation_df <- data.frame()
  for (j in 1:length(microbiome_mat)) {
    #j <- 1
    for (k in 1:length(Benthos_mat)) {
      #k <- 17    
      b = colnames(microbiome_mat)[j]
      c = colnames(Benthos_mat)[k]
      #str = paste0(a,"_",b,"_",c)/
      #print(str)
      mediation_tmp <- func_main(ALl_mat,a,b,c)
      mediation_df <- rbind(mediation_df,mediation_tmp)
    }
  }
  mediation_df
}

saveRDS(re_df,"~/Desktop/marine_sediment/results/Mediation_specific/pollution_microbiome_benthos1.rds")

## pollution -> function -> animal ---------------------------------------------
function_sign_df <- Pollution_path_list$COG_PATH_relist$sign_df %>%
  dplyr::rename(ID = feature, tmp = name) %>%
  left_join(Pollution_path_list$COG_PATH_list$Feature_ID)
function_mat <- Function_data$COG_list$COG_pathway_cpm %>%
  as.data.frame() %>%
  filter(Feature %in% function_sign_df$name) %>%
  column_to_rownames("Feature") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID")

all(rownames(Benthos_mat) == rownames(function_mat))

ALl_mat <- cbind(function_mat,Benthos_mat,Pollution_mat)

registerDoParallel(22)

re_df <- foreach (i = 1:dim(Pollution_mat)[2],.combine=rbind) %dopar% {
  #i <- 3
  set.seed(i)
  #print(i)
  a = colnames(Pollution_mat)[i]
  mediation_df <- data.frame()
  for (j in 1:length(function_mat)) {
    #j <- 1
    for (k in 1:length(Benthos_mat)) {
      #k <- 17    
      b = colnames(function_mat)[j]
      c = colnames(Benthos_mat)[k]
      #str = paste0(a,"_",b,"_",c)/
      #print(str)
      mediation_tmp <- func_main(ALl_mat,a,b,c)
      mediation_df <- rbind(mediation_df,mediation_tmp)
    }
  }
  mediation_df
}

saveRDS(re_df,"~/Desktop/marine_sediment/results/Mediation_specific/pollution_COG_benthos.rds")

## pollution -> microbiome -> animal -------------------------------------------
registerDoParallel(22)

re_df <- foreach (i = 1:dim(Pollution_mat)[2],.combine=rbind) %dopar% {
  #i <- 3
  set.seed(i)
  print(i)
  a = colnames(Pollution_mat)[i]
  mediation_df <- data.frame()
  for (j in 1:length(Benthos_mat)) {
    #j <- 1
    for (k in 1:length(microbiome_mat)) {
      #k <- 17    
      b = colnames(Benthos_mat)[j]
      c = colnames(microbiome_mat)[k]
      #str = paste0(a,"_",b,"_",c)/
      #print(str)
      mediation_tmp <- func_main(ALl_mat,a,b,c)
      mediation_df <- rbind(mediation_df,mediation_tmp)
    }
  }
  mediation_df
}

saveRDS(re_df,"~/Desktop/marine_sediment/results/Mediation_specific/pollution_benthos_microbiome.rds")
