## Load library   --------------------------------------------------------------
library(vegan)
library(ade4)
library(energy)
library(permute)
library(matrixStats)
source("/Users/linlin/software/mediation/MODIMA-master/modima.R")
library(tidyverse)
## OTU tables  -----------------------------------------------------------------
data_path <- "~/Desktop/Marine_Sediment/scripts_manuscript_new/shelby_data/mediation new/"
COG_pathway <- read_tsv(paste0(data_path,"COG_pathway_stat.flt.hellinger.tsv"))
COG <- read_tsv(paste0(data_path,"COG_cog.flt.hellinger.tsv"))
#Small_animal_Species_df <- read_tsv(paste0(data_path,"benthos_species.tsv"))
species_file <- paste0(data_path,"microbiome_species.otu")
read_taxon <- function(data){
  taxon <- read_tsv(data) %>%
    dplyr::select(.,-PC_S_51,-PC_W_33) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as.tibble() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.tibble() %>%
    dplyr::mutate_at(c(2:length(.)), as.numeric)
  ## 50% prevalence filter in at leat one site ---------------------------------
  Prev_df <- taxon %>%
    tidyr::separate(SampleID, into = c("Site", "Season","others"), sep = "_") %>%
    dplyr::select(.,-others,-Season) %>%
    tidyr::pivot_longer(2:(length(.)-1),names_to = "microbiome",values_to = "relativeAbundance") %>%
    dplyr::group_by(microbiome,Site) %>%
    dplyr::summarise_at("relativeAbundance", function(x) length(which(x!=0))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(microbiome) %>%
    dplyr::summarise_at("relativeAbundance", function(x) length(which(x>=3))) %>%
    dplyr::filter(relativeAbundance > 0)
  taxon <- taxon %>%
    tidyr::pivot_longer(2:length(.),names_to = "microbiome",values_to = "relativeAbundance") %>%
    dplyr::filter(microbiome %in% Prev_df$microbiome) %>%
    pivot_wider(names_from = microbiome,values_from = relativeAbundance) %>%
    column_to_rownames("SampleID")
  return(taxon)
}
Microbiome_species_otu <- read_taxon(species_file)
#small_animal_taxa <- read_csv(paste0(data_path,"TAXA_Linlin.csv"))
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
#Animal_TAXA <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/TAXA_Linlin.csv")
Animal_TAXA_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/AssignedTaxaTab.csv")
Animal_TAXA <- Animal_TAXA_file %>%
  filter(SourceID %in% "95_Match") %>%
  dplyr::select(.,-`...1`,-indtype,-SourceID) %>%
  rename(ID = Seq) %>%
  relocate(ID,.before = Kingdom)
#colnames(Animal_TAXA)[1] <- "ID"
dim(Animal_TAXA)
Animal_TAXA_df <- Animal_TAXA %>%
  unite(label,Kingdom:Species)
Animal_OTU_df <- Animal_OTU %>%
  left_join(Animal_SAMP) %>%
  relocate(c(Station_ID,Month_Ret), .after = ID) %>%
  mutate(Month_Ret = str_remove_all(Month_Ret,"[INUM]")) %>%
  mutate(ID = str_remove_all(ID,"ARMS_")) %>%
  unite("SampleID",c("Station_ID","Month_Ret","ID"))
Animal_otu <- Animal_OTU_df %>%
  #filter(SampleID %in% rownames(COG_path)) %>%
  column_to_rownames("SampleID") 
Small_animal_Species_df <- Animal_OTU_df %>%
  filter(SampleID %in% rownames(Microbiome_species_otu)) %>%
  arrange(SampleID) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(rowname %in% Animal_TAXA$ID) %>%
  # left_join(Animal_TAXA,by = c("rowname"= "ID")) %>%
  # filter(!is.na(Species)) %>%
  # group_by(Species) %>%
  # summarise_if(is.double,max) %>%
  # column_to_rownames("Species") %>%
  column_to_rownames() %>%
  t()  %>%
  as.data.frame() %>%
  rownames_to_column("ID")
# ## Significant features --------------------------------------------------------
# ancom2_gradient_list <- readRDS(paste0(data_path,"ancom2_gradient_list.rds"))
# 
# Microbiome_feature_tmp <- ancom2_gradient_list$Species$re_1$df %>%
#   dplyr::filter(!W %in% "Inf") %>%
#   dplyr::filter(detected_0.7 %in% "TRUE") %>% 
#   .$taxa_id
# Microbiome_feature <- ancom2_gradient_list$Species$re_1$res$eff_size %>% 
#   as.data.frame() %>%
#   dplyr::filter(. < 0) %>%
#   rownames_to_column() %>%
#   dplyr::filter(rowname %in% Microbiome_feature_tmp) %>%
#   .$rowname
# Microbiome_feature_pos_with_high_pollute <- setdiff(Microbiome_feature_tmp,Microbiome_feature)
# 
# Small_animal_org_reg_list <- readRDS(paste0(data_path,"ordinal_reg_list.rds"))
# #Animal_TAXA <- read_csv(paste0(data_path,"TAXA_Linlin.csv"))
# Animal_TAXA_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/AssignedTaxaTab.csv")
# Animal_TAXA <- Animal_TAXA_file %>%
#   filter(SourceID %in% "95_Match") %>%
#   dplyr::select(.,-`...1`,-indtype,-SourceID) %>%
#   rename(ID = Seq) %>%
#   relocate(ID,.before = Kingdom)
# #colnames(Animal_TAXA)[1] <- "ID"
# Animal_TAXA <- Animal_TAXA %>%
#   dplyr::select(.,-ID,-Species) %>%
#   dplyr::distinct()
# Small_animal_feature_genus <- Small_animal_org_reg_list$Genus$clm_df %>%
#   dplyr::filter(estimate > 0) %>%
#   dplyr::filter(name %in% Small_animal_org_reg_list$Genus$re_df_01$name) %>%
#   dplyr::select(name) %>%
#   dplyr::left_join(Animal_TAXA,by = c("name" = "Genus"))
# 
# Small_animal_feature_tmp <- Small_animal_org_reg_list$Species$clm_df %>%
#   dplyr::filter(estimate < 0) %>%
#   dplyr::filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
# Small_animal_feature <- Small_animal_feature_tmp$name
# Small_animal_feature_tmp_pos <- Small_animal_org_reg_list$Species$clm_df %>%
#   dplyr::filter(estimate > 0) %>%
#   dplyr::filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
# Small_animal_feature_pos <- Small_animal_feature_tmp_pos$name
# COG_list <- readRDS(paste0(data_path,"COG_ordinal_reg_list_1.rds"))
# COG_pathway_feature <- COG_list$COG_pathway_list$clm_df %>%
#   dplyr::filter(estimate < 0) %>%
#   dplyr::filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
# COG_pathway_feature_pos_with_high_pollute <- COG_list$COG_pathway_list$clm_df %>%
#   dplyr::filter(estimate > 0) %>%
#   dplyr::filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
## Load pollution data  ----------------------------------------------
pollution_file <- read_tsv(paste0(data_path,"year2019_pollution.tsv"))
pollution <- pollution_file %>%
  dplyr::relocate(c("Nitrate Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Total Inorganic Nitrogen (mg/L)","Total Nitrogen (mg/L)"),
           .after = "season") %>%
  dplyr::select(.,-year) %>%
  dplyr::rename(Site = StationID) %>%
  dplyr::rename(Season = season)
colnames(pollution) <- stringr::str_remove_all(colnames(pollution)," \\(μg/kg\\)") %>%
  stringr::str_remove_all(.," \\(mg/kg\\)") %>%
  stringr::str_remove_all(.," \\(mg/L\\)") %>%
  stringr::str_remove_all(.," \\(μg/L\\)")
Pollution_group <- read_tsv(paste0(data_path,"pollution_class.tsv"))
Pollution_group$name <- stringr::str_remove_all(Pollution_group$name," \\(μg/kg\\)") %>%
  stringr::str_remove_all(.," \\(mg/kg\\)") %>%
  stringr::str_remove_all(.," \\(mg/L\\)") %>%
  stringr::str_remove_all(.," \\(μg/L\\)")
normalize <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min)
  maxAttr=apply(x, 2, max)
  x <- sweep(x, 2, minAttr, FUN="-") 
  x=sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'normalized:min') = minAttr
  attr(x, 'normalized:max') = maxAttr
  return (x)
} 
pollution_mat <- data.frame("SampleID" = rownames(Microbiome_species_otu)) %>%
  dplyr::mutate(tmp = SampleID) %>%
  tidyr::separate(tmp, into = c("Site","Season","num")) %>%
  dplyr::mutate(Season = stringr::str_replace_all(Season,c("S"= "Summer","W" = "Winter"))) %>%
  dplyr::select(.,-num) %>%
  dplyr::left_join(pollution) %>%
  dplyr::arrange(SampleID) %>%
  column_to_rownames("SampleID") %>%
  dplyr::select(.,-Site,-Season) %>%
  normalize(.)
# ## Significant features matrix ----------------------------------------------------------
# microbiome_mat <- Microbiome_species_otu %>%
#   dplyr::select(Microbiome_feature_tmp) %>%
#   rownames_to_column("ID") %>%
#   dplyr::arrange(ID) %>%
#   column_to_rownames("ID")
# cog_path_mat <- COG_pathway %>%
#   dplyr::select(rowname,COG_pathway_feature$name,
#                 COG_pathway_feature_pos_with_high_pollute$name) %>%
#   dplyr::rename(ID = rowname) %>%
#   dplyr::arrange(ID) %>%
#   column_to_rownames("ID")
# small_animal_mat_neg <- Small_animal_Species_df %>%
#   dplyr::select(ID,Small_animal_feature) %>%
#   dplyr::arrange(ID) %>%
#   dplyr::filter(!ID %in% c("CI_S_35","PC_S_51","PC_W_33","SW_S_47")) %>%
#   column_to_rownames("ID")
# 
# small_animal_mat_pos <- Small_animal_Species_df %>%
#   dplyr::select(ID,Small_animal_feature_pos) %>%
#   dplyr::arrange(ID) %>%
#   filter(!ID %in% c("CI_S_35","PC_S_51","PC_W_33","SW_S_47")) %>%
#   column_to_rownames("ID")

## Mediation Analyses  pollution -> microbiome -> animal ----------------------------------------------------------

Bethos_mat <- Small_animal_Species_df %>%
  dplyr::filter(ID %in% rownames(Microbiome_species_otu)) %>%
  dplyr::arrange(ID) %>%
  column_to_rownames("ID") %>%
  as.matrix()

Microbiome_mat <- Microbiome_species_otu %>%
  rownames_to_column("ID") %>%
  dplyr::arrange(ID) %>%
  column_to_rownames("ID") %>%
  as.matrix()

COG_mat <- COG %>%
  dplyr::arrange(rowname) %>%
  column_to_rownames() %>%
  as.matrix()
 
## Mediation analyses  ---------------------------------------------------------
#exposure
dist_pollut = vegan::vegdist(pollution_mat, method="euclidean")
#response
dist_benthos = vegan::vegdist(Bethos_mat, method="jaccard")
#mediator
dist_taxon =vegan::vegdist(Microbiome_mat, method="bray")
dist_COG = vegan::vegdist(COG_mat, method="bray")

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
## dbRDA -----------------------------------------------------------------------
#load("~/Desktop/marine_sediment/RDS_files/beta_diversity_dbRDA_2019morevariables_list.RData")
RDA_pollut_feature_uni <- c("Arsenic","Copper","Zinc","Lead",
                        "Nitrate Nitrogen","Total Inorganic Nitrogen","Total Nitrogen",
                        "Benzo(ghi)perylene","Benzo(b)fluoranthene","Benzo(k)fluoranthene",
                        "Total Carbon (%w/w)")

## Distances   ---------------------------------------------------------
pollution_mat_dbRDA_uni <- pollution_mat %>%
  as.data.frame() %>%
  dplyr::select(RDA_pollut_feature_uni) %>%
  as.matrix()
#exposure
dist_pollut = vegan::vegdist(pollution_mat_dbRDA_uni, method="euclidean")
#response
dist_benthos = vegan::vegdist(Bethos_mat, method="jaccard")
#mediator
dist_taxon =vegan::vegdist(Microbiome_mat, method="bray")
dist_COG = vegan::vegdist(COG_mat, method="bray")
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

result_Pollution_Taxa_Benthos_list_dbRDA <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_Taxa_Benthos_dbRDA"
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

result_Pollution_COG_Benthos_list_dbRDA <- data.frame(
  Med_statistic = modima_results$statistic %>% as.numeric(),
  Med_percentage = (beta_MR_E * beta_EM)/beta_ER,
  Med_P = modima_results$p.value,
  beta_EM = beta_EM,
  P_EM = EM_list$p.value,
  beta_MR_E = beta_MR_E,
  P_MR_E = MR_E_list$p.value,
  beta_ER = beta_ER,
  P_ER = ER_list$p.value,
  name = "Pollution_COG_Benthos_dbRDA"
)
## Summary ---------------------------------------------------------------------
All_df <- rbind(
  result_Pollution_Taxa_Benthos_list_all,
  result_Pollution_COG_Benthos_list_all,
  result_Pollution_Taxa_Benthos_list_dbRDA,
  result_Pollution_COG_Benthos_list_dbRDA
)
write_tsv(All_df,"~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/mediation_new1/result_new_taxa_spe_panel.tsv")
