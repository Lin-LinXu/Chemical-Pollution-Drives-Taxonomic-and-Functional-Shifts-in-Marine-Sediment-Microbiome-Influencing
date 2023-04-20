## Load library ----------------------------------------------------------------
library(tidyverse)
library(mediation)
set.seed(123)
## OTU tables  -----------------------------------------------------------------
data_path <- "~/Projects/Marine_Sediment/mediation/mediation_new/"
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
    mutate_at(c(2:length(.)), as.numeric)
  ## 50% prevalence filter in at leat one site ---------------------------------
  Prev_df <- taxon %>%
    separate(SampleID, into = c("Site", "Season","others"), sep = "_") %>%
    dplyr::select(.,-others,-Season) %>%
    pivot_longer(2:(length(.)-1),names_to = "microbiome",values_to = "relativeAbundance") %>%
    group_by(microbiome,Site) %>%
    summarise_at("relativeAbundance", function(x) length(which(x!=0))) %>%
    ungroup() %>%
    group_by(microbiome) %>%
    summarise_at("relativeAbundance", function(x) length(which(x>=3))) %>%
    filter(relativeAbundance > 0)
  taxon <- taxon %>%
    pivot_longer(2:length(.),names_to = "microbiome",values_to = "relativeAbundance") %>%
    filter(microbiome %in% Prev_df$microbiome) %>%
    pivot_wider(names_from = microbiome,values_from = relativeAbundance) %>%
    column_to_rownames("SampleID")
  return(taxon)
}
Microbiome_species_otu <- read_taxon(species_file)
#small_animal_taxa <- read_csv(paste0(data_path,"TAXA_Linlin.csv"))
### Main part -----------------------------------------------------------------
Animal_OTU_file <- read_csv("~/Projects/Marine_Sediment/mediation/mediation_new/OTU_2019.csv")
colnames(Animal_OTU_file)[1] <- "ID"
Animal_OTU <- Animal_OTU_file %>%
  column_to_rownames("ID") %>%
  type_convert()
Animal_OTU[Animal_OTU > 0] <- 1
Animal_OTU <- Animal_OTU %>%
  rownames_to_column("ID")
Animal_SAMP <- read_csv("~/Projects/Marine_Sediment/mediation/mediation_new/SAMP_Linlin.csv")
colnames(Animal_SAMP)[1] <- "ID"
dim(Animal_SAMP)
Animal_SAMP <- Animal_SAMP %>%
  dplyr::select(.,-ARMS_ID,-Ret_ID)
Animal_TAXA_file <- read_csv("~/Projects/Marine_Sediment/mediation/mediation_new/Taxa_2019.csv")
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
  #filter(SampleID %in% rownames(COG_path)) %>%
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  left_join(Animal_TAXA,by = c("rowname"= "ID")) %>%
  filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarise_if(is.double,max) %>%
  column_to_rownames("Species") %>%
  t()  %>%
  as.data.frame() %>%
  rownames_to_column("ID")
## Significant features --------------------------------------------------------
ancom2_gradient_list <- readRDS(paste0(data_path,"ancom2_gradient_list.rds"))

Microbiome_feature_tmp <- ancom2_gradient_list$Species$re_1$df %>%
  filter(!W %in% "Inf") %>%
  filter(detected_0.7 %in% "TRUE") %>% 
  .$taxa_id
Microbiome_feature <- ancom2_gradient_list$Species$re_1$res$eff_size %>% 
  as.data.frame() %>%
  filter(. < 0) %>%
  rownames_to_column() %>%
  filter(rowname %in% Microbiome_feature_tmp) %>%
  .$rowname
Microbiome_feature_pos_with_high_pollute <- setdiff(Microbiome_feature_tmp,Microbiome_feature)

Small_animal_org_reg_list <- readRDS(paste0(data_path,"ordinal_reg_list.rds"))
Animal_TAXA <- read_csv(paste0(data_path,"Taxa_2019.csv"))
colnames(Animal_TAXA)[1] <- "ID"
Animal_TAXA <- Animal_TAXA %>%
  dplyr::select(.,-ID,-Species) %>%
  distinct()

Small_animal_feature_tmp <- Small_animal_org_reg_list$Species$clm_df %>%
  filter(estimate < 0) %>%
  filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
Small_animal_feature <- Small_animal_feature_tmp$name
Small_animal_feature_tmp_pos <- Small_animal_org_reg_list$Species$clm_df %>%
  filter(estimate > 0) %>%
  filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
Small_animal_feature_pos <- Small_animal_feature_tmp_pos$name
COG_list <- readRDS(paste0(data_path,"COG_ordinal_reg_list_1.rds"))
COG_pathway_feature <- COG_list$COG_pathway_list$clm_df %>%
  filter(estimate < 0) %>%
  filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
COG_pathway_feature_pos_with_high_pollute <- COG_list$COG_pathway_list$clm_df %>%
  filter(estimate > 0) %>%
  filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
## Load pollution data  ----------------------------------------------
pollution_file <- read_tsv(paste0(data_path,"year2019_pollution.tsv"))
pollution <- pollution_file %>%
  relocate(c("Nitrate Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Total Inorganic Nitrogen (mg/L)","Total Nitrogen (mg/L)"),
           .after = "season") %>%
  dplyr::select(.,-year) %>%
  rename(Site = StationID) %>%
  rename(Season = season)
colnames(pollution) <- str_remove_all(colnames(pollution)," \\(μg/kg\\)") %>%
  str_remove_all(.," \\(mg/kg\\)") %>%
  str_remove_all(.," \\(mg/L\\)") %>%
  str_remove_all(.," \\(μg/L\\)")
Pollution_group <- read_tsv(paste0(data_path,"pollution_class.tsv"))
Pollution_group$name <- str_remove_all(Pollution_group$name," \\(μg/kg\\)") %>%
  str_remove_all(.," \\(mg/kg\\)") %>%
  str_remove_all(.," \\(mg/L\\)") %>%
  str_remove_all(.," \\(μg/L\\)")
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
  mutate(tmp = SampleID) %>%
  separate(tmp, into = c("Site","Season","num")) %>%
  mutate(Season = str_replace_all(Season,c("S"= "Summer","W" = "Winter"))) %>%
  dplyr::select(.,-num) %>%
  left_join(pollution) %>%
  arrange(SampleID) %>%
  column_to_rownames("SampleID") %>%
  dplyr::select(.,-Site,-Season) %>%
  normalize(.)
## Mediation Analyses  pollution -> microbiome -> animal ----------------------------------------------------------
microbiome_mat <- Microbiome_species_otu %>%
  dplyr::select(Microbiome_feature_tmp) %>%
  rownames_to_column("ID") %>%
  arrange(ID) %>%
  column_to_rownames("ID")
cog_path_mat <- COG_pathway %>%
  dplyr::select(rowname,COG_pathway_feature$name,
                COG_pathway_feature_pos_with_high_pollute$name) %>%
  rename(ID = rowname) %>%
  arrange(ID) %>%
  column_to_rownames("ID")
small_animal_mat_neg <- Small_animal_Species_df %>%
  dplyr::select(ID,Small_animal_feature) %>%
  arrange(ID) %>%
  filter(!ID %in% c("CI_S_35","PC_S_51","PC_W_33","SW_S_47")) %>%
  column_to_rownames("ID")

small_animal_mat_pos <- Small_animal_Species_df %>%
  dplyr::select(ID,Small_animal_feature_pos) %>%
  arrange(ID) %>%
  filter(!ID %in% c("CI_S_35","PC_S_51","PC_W_33","SW_S_47")) %>%
  column_to_rownames("ID")

ALl_mat_neg <- cbind(microbiome_mat,cog_path_mat,small_animal_mat_neg,pollution_mat)
ALl_mat_pos <- cbind(microbiome_mat,cog_path_mat,small_animal_mat_pos,pollution_mat)

library(foreach)
library(doParallel)

## pollution -> microbiome -> animal (neg) -------------------------------------------
registerDoParallel(22)

re_df <- foreach (i = 1:dim(pollution_mat)[2],.combine=rbind) %dopar% {
  #i <- 3
  set.seed(i)
  #print(i)
  a = colnames(pollution_mat)[i]
  mediation_df <- data.frame()
  for (j in 1:length(microbiome_mat)) {
    #j <- 1
    for (k in 1:length(small_animal_mat_neg)) {
      #k <- 17
      b = colnames(microbiome_mat)[j]
      c = colnames(small_animal_mat_neg)[k]
      #str = paste0(a,"_",b,"_",c)
      #print(str)
      mediation_tmp <- func_main(ALl_mat_neg,a,b,c)
      mediation_df <- rbind(mediation_df,mediation_tmp)
    }
  }
  mediation_df
}

saveRDS(re_df,"~/Projects/Marine_Sediment/mediation/mediation_new/pollution_microbiome_benthos_neg.rds")
re_df <- readRDS("~/Projects/Marine_Sediment/mediation/mediation_new/pollution_microbiome_benthos_neg.rds")
re_df_tidy <- re_df %>%
  mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
  mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
  #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
  filter(!ACME_P_value %in% "") %>%
  mutate(BH_P = p.adjust(ACME_P_value)) %>%
  filter(ACME_P_value < 0.05,BH_P < 0.05) #%>%
  #filter(Beta_step1_iv < 0,Beta_step3_iv < 0)
  #filter(Beta_step2_iv > 0,Beta_step3_med > 0,Beta_step1_iv > 0,Beta_step3_iv >0)

registerDoParallel(22)

re_df_1 <- foreach (i = 1:dim(pollution_mat)[2],.combine=rbind) %dopar% {
  #i <- 3
  set.seed(i)
  #print(i)
  a = colnames(pollution_mat)[i]
  mediation_df <- data.frame()
  for (j in 1:length(small_animal_mat_neg)) {
    #j <- 1
    for (k in 1:length(microbiome_mat)) {
      #k <- 17
      b = colnames(small_animal_mat_neg)[j]
      c = colnames(microbiome_mat)[k]
      #str = paste0(a,"_",b,"_",c)
      #print(str)
      mediation_tmp <- func_main(ALl_mat_neg,a,b,c)
      mediation_df <- rbind(mediation_df,mediation_tmp)
    }
  }
  mediation_df
}

saveRDS(re_df_1,"~/Projects/Marine_Sediment/mediation/mediation_new/pollution_benthos_neg_microbiome.rds")
re_df_1 <- readRDS("~/Projects/Marine_Sediment/mediation/mediation_new/pollution_benthos_neg_microbiome.rds")
re_df_tidy1 <- re_df_1 %>%
  mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
  mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
  #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
  filter(!ACME_P_value %in% "") %>%
  mutate(BH_P = p.adjust(ACME_P_value)) %>%
  filter(ACME_P_value < 0.05,BH_P < 0.05)

re_df_tmp1 <- re_df_tidy1 %>%
  dplyr::select(iv,med,dv,ACME_P_value) %>%
  `colnames<-`(c("iv","dv","med","rev_p"))
re_df_tmp <- re_df_tidy %>%
  left_join(re_df_tmp1) %>%
  filter(!rev_p %in% 0)

write_tsv(re_df_tmp,"~/Projects/Marine_Sediment/mediation/mediation_new/mediation_species.tsv")
# re_df_tidy_pos_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv > 0)
# write_tsv(re_df_tidy_pos_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_neg_microbiome_pos.tsv")
# re_df_tidy_neg_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv < 0)
# write_tsv(re_df_tidy_neg_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_neg_microbiome_neg.tsv")


# ## pollution -> microbiome -> animal (pos) -------------------------------------------
# 
# registerDoParallel(22)
# 
# re_df <- foreach (i = 1:dim(pollution_mat)[2],.combine=rbind) %dopar% {
#   #i <- 3
#   set.seed(i)
#   #print(i)
#   a = colnames(pollution_mat)[i]
#   mediation_df <- data.frame()
#   for (j in 1:length(microbiome_mat)) {
#     #j <- 1
#     for (k in 1:length(small_animal_mat_pos)) {
#       #k <- 17
#       b = colnames(microbiome_mat)[j]
#       c = colnames(small_animal_mat_pos)[k]
#       #str = paste0(a,"_",b,"_",c)
#       #print(str)
#       mediation_tmp <- func_main(ALl_mat_pos,a,b,c)
#       mediation_df <- rbind(mediation_df,mediation_tmp)
#     }
#   }
#   mediation_df
# }
# 
# saveRDS(re_df,"~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_pos.rds")
# re_df <- readRDS("~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_pos.rds")
# 
# re_df_tidy <- re_df %>%
#   mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
#   mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
#   #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
#   filter(!ACME_P_value %in% "") %>%
#   mutate(BH_P = p.adjust(ACME_P_value)) %>%
#   filter(ACME_P_value < 0.05,BH_P < 0.05) %>%
#   filter(Beta_step1_iv > 0,Beta_step3_iv > 0)
# #filter(Beta_step2_iv > 0,Beta_step3_med > 0,Beta_step1_iv > 0,Beta_step3_iv >0)
# 
# re_df_tidy_pos_microbiome_pos_animal <- re_df_tidy %>%
#   filter(Beta_step2_iv > 0)
# write_tsv(re_df_tidy_pos_microbiome_pos_animal,"~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_pos_microbiome_pos.tsv")
# re_df_tidy_neg_microbiome_pos_animal <- re_df_tidy %>%
#   filter(Beta_step2_iv < 0)
# write_tsv(re_df_tidy_neg_microbiome_pos_animal,"~/Projects/Marine_Sediment/mediation/mediation/pollution_microbiome_benthos_pos_microbiome_neg.tsv")
# ## Mediation Analyses  pollution -> microbiome -> animal ----------------------------------------------------------
# microbiome_mat_neg <- Microbiome_species_otu %>%
#   dplyr::select(Microbiome_feature) %>%
#   rownames_to_column("ID") %>%
#   arrange(ID) %>%
#   column_to_rownames("ID")
# microbiome_mat_pos <- Microbiome_species_otu %>%
#   dplyr::select(Microbiome_feature_pos_with_high_pollute) %>%
#   rownames_to_column("ID") %>%
#   arrange(ID) %>%
#   column_to_rownames("ID")
# cog_path_mat <- COG_pathway %>%
#   dplyr::select(rowname,COG_pathway_feature$name,
#                 COG_pathway_feature_pos_with_high_pollute$name) %>%
#   rename(ID = rowname) %>%
#   arrange(ID) %>%
#   column_to_rownames("ID")
# small_animal_mat <- Small_animal_Species_df %>%
#   dplyr::select(ID,Small_animal_feature,Small_animal_feature_pos) %>%
#   arrange(ID) %>%
#   filter(!ID %in% c("CI_S_35","PC_S_51","PC_W_33","SW_S_47")) %>%
#   column_to_rownames("ID")
# 
# ALl_mat_neg <- cbind(microbiome_mat_neg,cog_path_mat,small_animal_mat,pollution_mat)
# ALl_mat_pos <- cbind(microbiome_mat_pos,cog_path_mat,small_animal_mat,pollution_mat)
# 
# library(foreach)
# library(doParallel)
# ## pollution -> animal -> microbiome (neg) -------------------------------------------
# registerDoParallel(22)
# 
# re_df <- foreach (i = 1:dim(pollution_mat)[2],.combine=rbind) %dopar% {
#   #i <- 3
#   set.seed(i)
#   #print(i)
#   a = colnames(pollution_mat)[i]
#   mediation_df <- data.frame()
#   for (j in 1:length(small_animal_mat)) {
#     #j <- 1
#     for (k in 1:length(microbiome_mat_neg)) {
#       #k <- 17
#       b = colnames(small_animal_mat)[j]
#       c = colnames(microbiome_mat_neg)[k]
#       #str = paste0(a,"_",b,"_",c)
#       #print(str)
#       mediation_tmp <- func_main(ALl_mat_neg,a,b,c)
#       mediation_df <- rbind(mediation_df,mediation_tmp)
#     }
#   }
#   mediation_df
# }
# 
# saveRDS(re_df,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_neg.rds")
# re_df <- readRDS("~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_neg.rds")
# re_df_tidy <- re_df %>%
#   mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
#   mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
#   #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
#   filter(!ACME_P_value %in% "") %>%
#   mutate(BH_P = p.adjust(ACME_P_value)) %>%
#   filter(ACME_P_value < 0.05,BH_P < 0.05) %>%
#   filter(Beta_step1_iv < 0,Beta_step3_iv < 0)
# #filter(Beta_step2_iv > 0,Beta_step3_med > 0,Beta_step1_iv > 0,Beta_step3_iv >0)
# 
# re_df_tidy_pos_animal_neg_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv > 0)
# write_tsv(re_df_tidy_pos_animal_neg_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_neg_benthos_pos.tsv")
# re_df_tidy_neg_animal_neg_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv < 0)
# write_tsv(re_df_tidy_neg_animal_neg_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_neg_benthos_neg.tsv")
# 
# ## pollution -> animal -> microbiome (pos) -------------------------------------------
# registerDoParallel(22)
# 
# re_df <- foreach (i = 1:dim(pollution_mat)[2],.combine=rbind) %dopar% {
#   #i <- 3
#   set.seed(i)
#   #print(i)
#   a = colnames(pollution_mat)[i]
#   mediation_df <- data.frame()
#   for (j in 1:length(small_animal_mat)) {
#     #j <- 1
#     for (k in 1:length(microbiome_mat_pos)) {
#       #k <- 17
#       b = colnames(small_animal_mat)[j]
#       c = colnames(microbiome_mat_pos)[k]
#       #str = paste0(a,"_",b,"_",c)
#       #print(str)
#       mediation_tmp <- func_main(ALl_mat_pos,a,b,c)
#       mediation_df <- rbind(mediation_df,mediation_tmp)
#     }
#   }
#   mediation_df
# }
# 
# saveRDS(re_df,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_pos.rds")
# re_df <- readRDS("~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_pos.rds")
# re_df_tidy <- re_df %>%
#   mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
#   mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
#   #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
#   filter(!ACME_P_value %in% "") %>%
#   mutate(BH_P = p.adjust(ACME_P_value)) %>%
#   filter(ACME_P_value < 0.05,BH_P < 0.05) %>%
#   filter(Beta_step1_iv > 0,Beta_step3_iv > 0)
# #filter(Beta_step2_iv > 0,Beta_step3_med > 0,Beta_step1_iv > 0,Beta_step3_iv >0)
# 
# re_df_tidy_pos_animal_pos_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv > 0)
# write_tsv(re_df_tidy_pos_animal_pos_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_pos_benthos_pos.tsv")
# re_df_tidy_neg_animal_pos_microbiome <- re_df_tidy %>%
#   filter(Beta_step2_iv < 0)
# write_tsv(re_df_tidy_neg_animal_pos_microbiome,"~/Projects/Marine_Sediment/mediation/mediation/pollution_benthos_microbiome_pos_benthos_neg.tsv")
# 

## Functions --------------------------------------------------------------------
func_main <- function(df,a,b,c)
{
  # df <-  ALl_mat
  # a <- "Benzo(a)anthracene"
  # b <- "Pirellula staleyi"
  # c <- "Syllidae_HK04"
  data_df <- func_df(df,a,b,c)
  func_step1_list <- func_step1(data_df)
  func_step2_list <- func_step2(data_df)
  if (func_step2_list$step2_re$`Pr(>|t|)` <= 0.05) {
    func_step3_list <- func_step3(data_df)
    Step3_P <- func_step3_list$step3_re %>% filter(terms %in% "med") %>% .$`Pr(>|t|)`
    Step3_beta <- func_step3_list$step3_re %>% filter(terms %in% "iv") %>% .$Estimate %>% abs()
    step1_beta <- func_step1_list$step1_re$Estimate %>% abs()
    Delta_beta <- step1_beta - Step3_beta
    if (Step3_P <= 0.05 && Delta_beta > 0) {
      func_step4_df <- func_step4(func_step2_list$fit.mediator,func_step3_list$fit.dv)
      re_df <- data.frame(
        iv = a,
        med = b,
        dv = c,
        Beta_step1_iv = func_step1_list$step1_re$Estimate,
        Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
        Beta_step2_iv = func_step2_list$step2_re$Estimate,
        Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
        Beta_step3_iv = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$Estimate,
        Beta_step3_iv_p = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
        Beta_step3_med = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$Estimate,
        Beta_step3_med_p = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$`Pr(>|t|)`,
        ACME = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$Estimate,
        ADE = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ADE") %>%
          .$Estimate,
        ACME_lower_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`95% CI Lower`,
        ACME_upper_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`95% CI Upper`,
        ACME_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ACME") %>%
          .$`p-value`,
        ADE_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "ADE") %>%
          .$`p-value`,
        Percentage_Medi_P_value = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`p-value`,
        Percentage_Medi = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$Estimate,
        Percentage_Medi_lower_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`95% CI Lower`,
        Percentage_Medi_upper_95_CI = func_step4_df %>% 
          as.data.frame() %>%
          rownames_to_column() %>% 
          filter(rowname %in% "Prop. Mediated") %>%
          .$`95% CI Upper`
      )
    }else
      {
      re_df <- data.frame(
        iv = a,
        med = b,
        dv = c,
        Beta_step1_iv = func_step1_list$step1_re$Estimate,
        Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
        Beta_step2_iv = func_step2_list$step2_re$Estimate,
        Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
        Beta_step3_iv = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$Estimate,
        Beta_step3_iv_p = func_step3_list$step3_re %>%
          filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
        Beta_step3_med = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$Estimate,
        Beta_step3_med_p = func_step3_list$step3_re %>%
          filter(terms %in% "med") %>% .$`Pr(>|t|)`,
        ACME = "",
        ADE = "",
        ACME_lower_95_CI = "",
        ACME_upper_95_CI = "",
        ACME_P_value = "",
        ADE_P_value = "",
        Percentage_Medi_P_value = "",
        Percentage_Medi = "",
        Percentage_Medi_lower_95_CI = "",
        Percentage_Medi_upper_95_CI = ""
      )
    }
  }else{
    func_step3_list <- func_step3(data_df)
    re_df <- data.frame(
      iv = a,
      med = b,
      dv = c,
      Beta_step1_iv = func_step1_list$step1_re$Estimate,
      Beta_step1_iv_p = func_step1_list$step1_re$`Pr(>|t|)`,
      Beta_step2_iv = func_step2_list$step2_re$Estimate,
      Beta_step2_iv_p = func_step2_list$step2_re$`Pr(>|t|)`,
      Beta_step3_iv = func_step3_list$step3_re %>%
        filter(terms %in% "iv") %>% .$Estimate,
      Beta_step3_iv_p = func_step3_list$step3_re %>%
        filter(terms %in% "iv") %>% .$`Pr(>|t|)`,
      Beta_step3_med = func_step3_list$step3_re %>%
        filter(terms %in% "med") %>% .$Estimate,
      Beta_step3_med_p = func_step3_list$step3_re %>%
        filter(terms %in% "med") %>% .$`Pr(>|t|)`,
      ACME = "",
      ADE = "",
      ACME_lower_95_CI = "",
      ACME_upper_95_CI = "",
      ACME_P_value = "",
      ADE_P_value = "",
      Percentage_Medi_P_value = "",
      Percentage_Medi = "",
      Percentage_Medi_lower_95_CI = "",
      Percentage_Medi_upper_95_CI = ""
    )
  }
  return(re_df)
}

func_df <- function(df,a,b,c)
{
 # df <-  ALl_mat
 # a <- "Nitrate Nitrogen"
 # b <- "Nitrospina gracilis"
 # c <- "Polynoidae HK03"
 data_df <- df %>%
   dplyr::select(a,b,c) %>%
   `colnames<-`(c("iv","med","dv"))
 return(data_df)
}

func_step1 <- function(df)
{
  #Step 1: The total effect
  fit.totaleffect=lm(dv~iv,df)
  #summary(fit.totaleffect)
  step1_re <- as.data.frame(coef(summary(fit.totaleffect))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% "iv") %>%
    mutate(model = "step1")
  return_list <- list(
    fit.totaleffect = fit.totaleffect,
    step1_re = step1_re
  )
  return(return_list)
}
  
func_step2 <- function(df)
{
  #Step 2: The effect of the IV onto the mediator
  fit.mediator=lm(med~iv,df)
  #summary(fit.mediator)
  step2_re <- as.data.frame(coef(summary(fit.mediator))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% "iv") %>%
    mutate(model = "step2")
  return_list <- list(
    fit.mediator = fit.mediator,
    step2_re = step2_re
  )
  return(return_list)
}

func_step3 <- function(df)
{
  #Step 3: The effect of the mediator on the dependent variable
  fit.dv=lm(dv~iv+med,df)
  #summary(fit.dv)
  step3_re <- as.data.frame(coef(summary(fit.dv))) %>%
    rownames_to_column("terms") %>%
    filter(terms %in% c("iv","med")) %>%
    mutate(model = "step3")
  return_list <- list(
    fit.dv = fit.dv,
    step3_re = step3_re
  )
  return(return_list)
}

func_step4 <- function(fit.mediator,fit.dv)
{
  #Step 4: Causal Mediation Analysis
  results = mediate(fit.mediator, fit.dv, treat='iv', mediator='med', boot=T)
  #summary(results)
  summary_tb <- extract_mediation_summary(summary(results))
  return(summary_tb)
}

## Sub functions ---------------------------------------------------------------
extract_mediation_summary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}

