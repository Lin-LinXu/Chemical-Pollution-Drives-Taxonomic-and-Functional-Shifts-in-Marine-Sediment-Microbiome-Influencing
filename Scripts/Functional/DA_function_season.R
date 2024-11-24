# Title: Differential Abundance microbiome function (season)
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

library(Maaslin2)
## Load data -------------------------------------------------------------------
data_PATH <- readRDS("/data/Desktop/marine_sediment/results/Function_DA_PATHWAY_pollution_List.rds")
data <- readRDS("/data/Desktop/marine_sediment/results/Function_DA_pollution_List.rds")
## function --------------------------------------------------------------------
func_Maaslin2_season <- function(data_list){
  fit_data = Maaslin2(input_data     = data_list$taxon, 
                      input_metadata = data_list$meta_data, 
                      min_prevalence = 0,
                      min_abundance = 0,
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOG",
                      output         = "demo_output1", 
                      fixed_effects  = c("Season"),
                      random_effects = c("Site")
  )
  
  directoryPath <- "/data/demo_output1"
  unlink(directoryPath, recursive = TRUE)
  
  library(tidyverse)
  
  sign_df <- fit_data$results %>%
    filter(metadata %in% "Season") %>%
    filter(pval <= 0.05, qval <= 0.20) %>%
    mutate()
  
  re_list <- list(
    results = fit_data$results,
    sign_df = sign_df
  )
  return(re_list)
}
## results ---------------------------------------------------------------------
COG_PATH_relist <- func_Maaslin2_season(data_PATH$COG_PATH_list)
COG_category_relist <- func_Maaslin2_season(data_PATH$COG_category_list)

COG_relist <- func_Maaslin2_season(data$COG_list)

Functional_DA_season_List <- list(
  COG_PATH_relist = COG_PATH_relist,
  COG_category_relist = COG_category_relist,
  COG_relist = COG_relist
)

saveRDS(Functional_DA_season_List,"/data/Desktop/marine_sediment/results/Functional_DA_season_List.rds")
