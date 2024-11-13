library(Maaslin2)
## Load data -------------------------------------------------------------------
data <- readRDS("/data/Desktop/marine_sediment/results/Benthos/Animal_List.rds")
## function --------------------------------------------------------------------
func_Maaslin2_season <- function(data_list){
  fit_data = Maaslin2(input_data     = data_list$taxon, 
                      input_metadata = data_list$meta_data, 
                      min_prevalence = 0,
                      min_abundance = 0,
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "NONE",
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
Animal_season_relist <- func_Maaslin2_season(data$Animal_otu_list)

Animal_List <- list(
  Animal_otu_relist = Animal_season_relist,
  Animal_otu_list = data$Animal_otu_list
)

Animal_otu_df <- Animal_season_relist$sign_df %>%
  left_join(data$Animal_otu_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

saveRDS(Animal_List,"/data/Desktop/marine_sediment/results/Benthos/Animal_seasonList.rds")

Animal_List <- readRDS("~/Desktop/marine_sediment/results/Benthos/Animal_seasonList.rds")

Animal_otu_df <- Animal_season_relist$sign_df %>%
  left_join(data$Animal_otu_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))
