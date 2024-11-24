# Title: Mediation LDM
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MI

## Load library ----------------------------------------------------------------
library(LDM)
library(tidyverse)
library(phyloseq)
## Load data  ------------------------------------------------------------------
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
species_df <- read_taxon(paste0(data_path,"species.otu"))
### Function data --------------------------------------------------------------
Functional_list <- readRDS("~/Desktop/marine_sediment/results/Functional_beta_diversity.rds")
### Pollution ------------------------------------------------------------------
file1 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ms5.csv")
file2 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ns6.csv")
file3 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ps6.csv")
file4 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ss1-ss3-ss5.csv")
file5 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ts2.csv")
file6 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5yrs-ms8.csv")
file_sediment <- rbind(file1,file2,file3,file4,file5,file6)
pollution_sediment <- read_csv2("~/Desktop/marine_sediment/data/pollution.csv") %>%
  select(STATION,StationID) %>%
  distinct() %>%
  rename(Station = STATION)
df_sediment <- file_sediment %>%
  #filter(!Station %in% "MS8") %>%
  left_join(pollution_sediment) %>%
  relocate(StationID,.before = Station)
df_tidy_sediment <- df_sediment %>%
  mutate_all(funs(str_replace(., "<", ""))) %>%
  mutate_all(funs(str_replace(., "N/A", NA_character_))) %>%
  as_tibble() %>%
  type_convert() %>%
  separate(Dates,c("year","month","day"), sep = "-") %>%
  mutate(month = str_replace(month, "0", "")) %>%
  mutate(
    season = case_when(
      month %in% 10:12 ~ "Fall",
      month %in%  1:3  ~ "Winter",
      month %in%  4:6  ~ "Spring",
      TRUE ~ "Summer")) %>%
  relocate(season,.before = month)

sediment_all_names <- read_tsv("~/Desktop/marine_sediment/renew_pollutionbin/sediment_pollute.tsv")
remove_pollutants <- sediment_all_names %>%
  filter(categories %in% c("Organic pollutants",
                           "physical indicators",
                           "Others")
  )
Score_sediment <- df_tidy_sediment %>%
  filter(year %in% c("2019"),season %in% c("Summer","Winter")) %>%
  type_convert() %>%
  group_by(StationID,season) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  dplyr::select(.,-c("year","month","Sample No")) %>%
  #column_to_rownames("StationID") %>%
  dplyr::rename(`TotalPCBs (μg/kg)` = `Total Polychlorinated Biphenyls (μg/kg)`) %>%
  dplyr::select(.,-c(`TotalPCBs (μg/kg)`,`Silver (mg/kg)`)) %>%
  dplyr::select(.,-remove_pollutants$name) %>%
  unite("group", StationID:season)
df_form <- data.frame(sample = rownames(species_df)) %>%
  mutate(tmp = sample) %>%
  separate(tmp, c("Site", "Season","others")) %>%
  mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
  dplyr::select(.,-others) %>%
  unite("group", Site:Season)
merge_tab <- df_form %>%
  left_join(Score_sediment)
rownames(merge_tab) <- merge_tab$sample
Pollution_df = merge_tab[,!(names(merge_tab) %in% c("group","sample","Site","Season"))] 
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

## Prepare dataframe -----------------------------------------------------------
Input_data <- list(
  Pollution_df = Pollution_df,
  Animal_otu = Animal_otu,
  species_df = species_df
)

saveRDS(Input_data,"~/Desktop/marine_sediment/scripts/Mediation/Mediation_panel/LDM-med/Input_data.rds")
exposure_ID <- data.frame(name = colnames(exposure_df)) %>%
  mutate(ID = paste0("exposure",seq(1,dim(exposure_df)[2])))

outcome_ID <- data.frame(name = colnames(outcome_df)) %>%
  mutate(ID = paste0("outcome",seq(1,dim(outcome_df)[2])))

Mediator_ID <- data.frame(name = colnames(Mediator_df)) %>%
  mutate(ID = paste0("Mediator",seq(1,dim(Mediator_df)[2])))

exposure_df <- Pollution_df %>% `colnames<-`(exposure_ID$ID)
outcome_df <- Animal_otu %>% `colnames<-`(outcome_ID$ID)
Mediator_df <- species_df %>% `colnames<-`(Mediator_ID$ID)



all(rownames(exposure_df) == rownames(outcome_df))
exposure_string <- paste(colnames(exposure_df), collapse = " + ")
outcome_string <- paste(colnames(outcome_df), collapse = " + ")

meta_df <- cbind(exposure_df,outcome_df)
## Mediation test --------------------------------------------------------------
formula_str <- as.formula(paste0("Mediator_df ~"," (",exposure_string, ") + ",
                             "(",outcome_string,")"))

test_str <- as.formula(paste0("Mediator_df ~"," (","exposure1 + exposure2", ") + ",
                              "(","outcome1 + outcome2",")"))

res.ldm.med <- ldm(formula=formula_str,
                     data=meta_df, seed=67817, n.cores=2,
                     test.mediation=TRUE)   # parameter for requesting LDM-omni3

res.ldm.med$med.p.global.omni 
res.ldm.med$med.detected.otu.omni 
