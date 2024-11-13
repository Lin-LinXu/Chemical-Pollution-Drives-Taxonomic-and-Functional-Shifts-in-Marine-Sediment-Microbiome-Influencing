## Load library   --------------------------------------------------------------
library(vegan)
library(ade4)
library(energy)
library(permute)
library(matrixStats)
library(phyloseq)
source("/Users/linlin/software/mediation/MODIMA-master/modima.R")
library(tidyverse)
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
genus_df <- read_taxon(paste0(data_path,"genus.otu")) 
species_df <- read_taxon(paste0(data_path,"species.otu"))
### Function data --------------------------------------------------------------
Functional_list <- readRDS("~/Desktop/marine_sediment/results/Functional_beta_diversity.rds")
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
## Distance matrix  ----------------------------------------------------------------
library(vegan)
spe.dist <- vegdist(species_df, method = 'bray')
genus.dist <- vegdist(genus_df, method = 'bray')
animal.dist <- vegdist(Animal_otu, method = 'jaccard')
cog.dist <- Functional_list$LogEuclidean_dist[[1]]
ko.dist <- Functional_list$LogEuclidean_dist[[2]]
pollution.dist <- vegdist(Pollution_df, method = 'euclidean')

all(colnames(as.matrix(spe.dist)) == colnames(as.matrix(pollution.dist)))
## Mediation Analyses  pollution -> microbiome -> animal ----------------------------------------------------------

#exposure
dist_pollut = pollution.dist
#response
dist_benthos = animal.dist
#mediator
dist_taxon = spe.dist
dist_COG = cog.dist

## Pollution -> Microbiome -> Benthos ------------------------------------------
set.seed(12345)
modima_results = modima(exposure=dist_pollut, 
                        mediator=dist_taxon, 
                        response=dist_benthos, 
                        nrep=9999)

set.seed(12345)
ER_M_list <- pdcor.test(dist_pollut, dist_benthos, dist_taxon, R=9999)
set.seed(12345)
MR_E_list <- pdcor.test(dist_taxon, dist_benthos, dist_pollut, R=9999)
set.seed(12345)
beta_MR_E <- pdcor(dist_taxon, dist_benthos, dist_pollut)  %>% as.numeric()
set.seed(12345)
beta_ER <- bcdcor(dist_pollut,dist_benthos)  %>% as.numeric()
set.seed(12345)
ER_list <- dcor.test(dist_pollut,dist_benthos,R = 9999)
set.seed(12345)
beta_EM <- bcdcor(dist_pollut,dist_taxon)  %>% as.numeric()
set.seed(12345)
EM_list <- dcor.test(dist_pollut,dist_taxon,R = 9999)

result_Pollution_Taxa_Benthos_list_all <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_Taxa_Benthos_all"
)

## Pollution -> Microbiome COG -> Benthos ------------------------------------------
set.seed(12345)
modima_results = modima(exposure=dist_pollut, 
                        mediator=dist_COG, 
                        response=dist_benthos, 
                        nrep=9999)

set.seed(12345)
ER_M_list <- pdcor.test(dist_pollut, dist_benthos, dist_COG, R=9999)
set.seed(12345)
MR_E_list <- pdcor.test(dist_COG, dist_benthos, dist_pollut, R=9999)
set.seed(12345)
beta_MR_E <- pdcor(dist_COG, dist_benthos, dist_pollut)  %>% as.numeric()
set.seed(12345)
beta_ER <- bcdcor(dist_pollut,dist_benthos)  %>% as.numeric()
set.seed(12345)
ER_list <- dcor.test(dist_pollut,dist_benthos,R = 9999)
set.seed(12345)
beta_EM <- bcdcor(dist_pollut,dist_COG)  %>% as.numeric()
set.seed(12345)
EM_list <- dcor.test(dist_pollut,dist_COG,R = 9999)

result_Pollution_COG_Benthos_list_all <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_COG_Benthos_all"
)

All_df <- rbind(
  result_Pollution_Taxa_Benthos_list_all,
  result_Pollution_COG_Benthos_list_all
)
write_tsv(All_df,"~/Desktop/marine_sediment/results/Benthos/mediation_panel/Panel_modima.tsv")

## Pollution -> Benthos -> Microbiome ------------------------------------------
set.seed(12345)
modima_results = modima(exposure=dist_pollut, 
                        mediator=dist_benthos, 
                        response=dist_taxon, 
                        nrep=9999)

set.seed(12345)
ER_M_list <- pdcor.test(dist_pollut, dist_taxon, dist_benthos, R=9999)
set.seed(12345)
MR_E_list <- pdcor.test(dist_benthos, dist_taxon, dist_pollut, R=9999)
set.seed(12345)
beta_MR_E <- pdcor(dist_benthos, dist_taxon, dist_pollut)  %>% as.numeric()
set.seed(12345)
beta_ER <- bcdcor(dist_pollut,dist_taxon)  %>% as.numeric()
set.seed(12345)
ER_list <- dcor.test(dist_pollut,dist_taxon,R = 9999)
set.seed(12345)
beta_EM <- bcdcor(dist_pollut,dist_benthos)  %>% as.numeric()
set.seed(12345)
EM_list <- dcor.test(dist_pollut,dist_benthos,R = 9999)

result_Pollution_Benthos_Taxa_list_all <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_Benthos_Taxa_all"
)

## Pollution -> Benthos --> Microbiome COG------------------------------------------
set.seed(12345)
modima_results = modima(exposure=dist_pollut, 
                        mediator=dist_benthos, 
                        response=dist_COG, 
                        nrep=9999)

set.seed(12345)
ER_M_list <- pdcor.test(dist_pollut, dist_COG, dist_benthos, R=9999)
set.seed(12345)
MR_E_list <- pdcor.test(dist_benthos, dist_COG, dist_pollut, R=9999)
set.seed(12345)
beta_MR_E <- pdcor(dist_benthos, dist_COG, dist_pollut)  %>% as.numeric()
set.seed(12345)
beta_ER <- bcdcor(dist_pollut,dist_COG)  %>% as.numeric()
set.seed(12345)
ER_list <- dcor.test(dist_pollut,dist_COG,R = 9999)
set.seed(12345)
beta_EM <- bcdcor(dist_pollut,dist_benthos)  %>% as.numeric()
set.seed(12345)
EM_list <- dcor.test(dist_pollut,dist_benthos,R = 9999)

result_Pollution_Benthos_COG_list_all <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_Benthos_COG_all"
)

All_df <- rbind(
  result_Pollution_Benthos_Taxa_list_all,
  result_Pollution_Benthos_COG_list_all
)
write_tsv(All_df,"~/Desktop/marine_sediment/results/Benthos/mediation_panel/Panel_modima_benthos.tsv")
