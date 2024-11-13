## Load library ----------------------------------------------------------------
library(tidyverse)
library(rio)
## Read files ------------------------------------------------------------------
Metal_list <- import_list("~/Desktop/marine_sediment/literature/new_function_articles/metal_resistence_KOs.xlsx")
Others_list <- import_list("~/Desktop/marine_sediment/literature/new_function_articles/Other_functions.xlsx")
## Tidy data -------------------------------------------------------------------
Metal_df <- Metal_list$`Table S1`
Others_df <- Others_list$Sheet1 %>%
  dplyr::rename(KO = `KEGG KO ID`,
                Functional_module = `Functional module`,
                Name = Annotation
                ) %>%
  filter(!Functional_module %in% c("Iron-related metabolism",
                                   "Manganese-related metabolism",
                                   "Flagellar assembly")) %>%
  dplyr::select(Functional_module,KO,Name)

#Metal_table <- table(Metal_df$Metal) %>% as.data.frame()

Metal_tidy_df <- Metal_df %>%
  separate_rows(Metal, sep = "/") %>%
  filter(!Metal %in% c("Ag", "Au", "Bi","Co","Mg","W","Mo","Te")) %>%
  mutate(Metal = if_else(Mechanisms == "Extracellular sequestration",
                         "Extracellular sequestration",Metal)) %>%
  filter(!Metal %in% c("-")) %>%
  dplyr::rename(Functional_module = Metal) %>%
  dplyr::select(Functional_module,KO,Name)

Final_KO_df <- list(
  Metal_df = Metal_tidy_df,
  Others_df = Others_df
)
saveRDS(Final_KO_df,"~/Desktop/marine_sediment/results/KO_pathway_files/Function_KO_groups.rds")

## Loading data ----------------------------------------------------------------
data_path <- "~/Desktop/marine_sediment/results/"
function_list <- readRDS(paste0(data_path,"functional_profile.rds"))
KO_df <- function_list$KO_list$KO_cpm
Meta_resistence_KO_cpm <- KO_df %>%
  filter(Feature %in% Final_KO_df$Metal_df$KO) %>%
  dplyr::rename(KO = Feature) %>%
  left_join(Final_KO_df$Metal_df) %>%
  dplyr::select(.,-c("KO","Name")) %>%
  ungroup() %>%
  group_by(Functional_module) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  dplyr::rename(Feature = Functional_module)
Others_KO_cpm <- KO_df %>%
  filter(Feature %in% Final_KO_df$Others_df$KO) %>%
  dplyr::rename(KO = Feature) %>%
  left_join(Final_KO_df$Others_df) %>%
  dplyr::select(.,-c("KO","Name")) %>%
  ungroup() %>%
  group_by(Functional_module) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  dplyr::rename(Feature = Functional_module)
KO_detail_function_list <- list(
  Meta_resistence_KO_cpm = Meta_resistence_KO_cpm,
  Others_KO_cpm = Others_KO_cpm
)

saveRDS(KO_detail_function_list,"~/Desktop/marine_sediment/results/KO_pathway_files/KO_detail_function_data.rds")

