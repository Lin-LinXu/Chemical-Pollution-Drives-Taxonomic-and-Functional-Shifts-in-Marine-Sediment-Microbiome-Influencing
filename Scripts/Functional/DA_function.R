## Loading required packages ---------------------------------------------------
library(Maaslin2)
library(tidyverse)
library(phyloseq)
#library(MatrixExtra)
## Loading data ----------------------------------------------------------------
data_path <- "/data/Desktop/marine_sediment/results/"
function_list <- readRDS(paste0(data_path,"functional_profile.rds"))
## Prepare df ------------------------------------------------------------------
Pre_filter <- function(data)
{
  #data <- function_list$KO_list$KO_pathway_cpm
  
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
COG_PATH_list <- Pre_filter(function_list$COG_list$COG_pathway_cpm)
all(rownames(COG_PATH_list$taxon) == rownames(COG_PATH_list$meta_data))

COG_list <- Pre_filter(function_list$COG_list$COG_cpm)
all(rownames(COG_list$taxon) == rownames(COG_list$meta_data))


COG_category_list <- Pre_filter(function_list$COG_list$COG_category_cpm)
all(rownames(COG_category_list$taxon) == rownames(COG_category_list$meta_data))
## Maaslin2 --------------------------------------------------------------------
func_Maaslin2 <- function(data_list){
  fit_data = Maaslin2(input_data     = data_list$taxon, 
                      input_metadata = data_list$meta_data, 
                      min_prevalence = 0,
                      min_abundance = 0,
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOG",
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
COG_PATH_relist <- func_Maaslin2(COG_PATH_list)
COG_category_relist <- func_Maaslin2(COG_category_list)

COG_PATH_df <- COG_PATH_relist$sign_df %>%
  left_join(COG_PATH_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))

COG_category_df <- COG_category_relist$sign_df %>%
  left_join(COG_category_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))


Function_PATH_List <- list(
  COG_PATH_relist = COG_PATH_relist,
  COG_category_relist = COG_category_relist,
  COG_PATH_list = COG_PATH_list,
  COG_category_list = COG_category_list
)

saveRDS(Function_PATH_List,"/data/Desktop/marine_sediment/results/Function_DA_PATHWAY_pollution_List.rds")
Function_PATH_List <- readRDS("~/Desktop/marine_sediment/results/Function_DA_PATHWAY_pollution_List.rds")
COG_relist <- func_Maaslin2(COG_list)

COG_df <- COG_relist$sign_df %>%
  left_join(COG_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))


Function_List <- list(
  COG_relist = COG_relist,
  COG_list = COG_list
)

saveRDS(Function_List,"/data/Desktop/marine_sediment/results/Function_DA_pollution_List.rds")
Function_List <- readRDS("~/Desktop/marine_sediment/results/Function_DA_pollution_List.rds")

  
COG_df <- Function_List$COG_relist$sign_df %>%
  left_join(Function_List$COG_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name)) %>%
  mutate(sign = if_else(coef > 0, "High","Low"))


library(xlsx)
write.xlsx(COG_df, "~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiom_resubmit_Cell_Reports/Supplementary_tables_new/Supplementary_Tabs.xlsx",
           sheetName="Table S6",row.names = F,append=TRUE)

