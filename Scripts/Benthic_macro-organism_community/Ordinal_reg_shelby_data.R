## Load library ----------------------------------------------------------------
library(ordinal)
library(tidyverse)
library(broom)
library(ComplexHeatmap)
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
Animal_TAXA_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/Taxa_2019.csv")
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
Animal_species <- Animal_OTU_df %>%
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
### Setting color parameters ---------------------------------------------------
cbPalette1 <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
envPalette1 <- c("#0072B2","#E69F00","#D55E00")
envPalette <- c("High" ="#D55E00","Medium" = "#E69F00","Low" = "#0072B2")
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
seasonPalette <- c("Summer" = "#0073C2FF","Winter" = "#EFC000FF")
## Function --------------------------------------------------------------------
# prop odds test
# see https://stat.ethz.ch/pipermail/r-help/2014-November/423706.html
test_prop_odds <- function(df) {
  #
  # test if the proportional odds assumption holds
  # 1. create a clm with equal slope and another clm with multiple slopes
  # 2. anova to test whether a difference exist
  # Proportional odds assumption holds if p > 0.05
  #
  # @param df data frame
  # @return Chi-Squared p-value of anova
  #
  p.val <- tryCatch({
    clm_prop_odds <- clm(Gradient ~ value + Season, link="logit", data=df)
    clm_multi_slope <- clm(Gradient ~ Season, link="logit",  nominal = ~value, data=df)
    anova(clm_prop_odds, clm_multi_slope)$`Pr(>Chisq)`[2]
    #}, warning = function(war) {
    #  0
  }, error = function(err) {
    0
  })
  return(as.data.frame(p.val))
}
Ordinal_reg <- function(df,vec){
  # df <- Genus_df
  # vec <- c("Low","Medium","High")
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  df <- df %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble() %>%
    left_join(gradient) %>%
    relocate(c("Site","Season","Gradient"),.after = ID) %>%
    pivot_longer(5:length(.)) %>%
    mutate(Gradient = factor(Gradient,levels = c("Low","Medium","High")))
  gradient_vec <- c(paste(vec[1],vec[2],sep = "|"),paste(vec[2],vec[3],sep = "|"))
  gradient_vec1 <- c(paste(gradient_vec[1],".value",sep = ""),
                     paste(gradient_vec[2],".value",sep = "") )
  
  clm_df <- df %>%
    ungroup() %>%
    group_by(name) %>% # test for each metabolite separatly
    #do(tidy(polr(Class2 ~ abundance + Age, data = ., Hess=TRUE))) %>%
    do(tidy(clm(Gradient ~ value + Season, link="logit", data = .))) %>% # ordinal regression
    ungroup() %>%
    group_by(term) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>% # p value correction
    ungroup() %>%
    filter(term %in% "value")
  
  significant_clm_cags <- clm_df %>% 
    filter(p.value < 0.05,q.value < 0.05 ) %>% 
    pluck("name")
  significant_clm_cags_01 <- clm_df %>% 
    filter(p.value < 0.05,q.value < 0.1 ) %>% 
    pluck("name")
  significant_clm_cags_025 <- clm_df %>% 
    filter(p.value < 0.05,q.value < 0.25 ) %>% 
    pluck("name")
  clm_df <- mutate(clm_df, is_clm_significant_cag = name %in% significant_clm_cags)
  # Test prop odds
  prop_odds_passed_cags <- df %>%
    ungroup() %>%
    group_by(name) %>%
    do(test_prop_odds(.)) %>%
    filter(p.val > 0.05) %>% # H0: no diff btw single and multi slope model
    pluck("name")
  
  significant_clm_cags_025_tmp <- significant_clm_cags_025
  signif_assoc_cags <- intersect(prop_odds_passed_cags, significant_clm_cags)
  signif_assoc_cags_01 <- intersect(prop_odds_passed_cags, significant_clm_cags_01)
  signif_assoc_cags_025 <- intersect(prop_odds_passed_cags, significant_clm_cags_025)
  left_over <- setdiff(significant_clm_cags_025_tmp,signif_assoc_cags_025)
  clm_df_unpropotional <- df %>%
    ungroup() %>%
    #filter(name %in% left_over) %>%
    group_by(name) %>% # test for each metabolite separatly
    summarise(tidy(clm(Gradient ~ Season,nominal = ~value,data = .))) %>% # ordinal regression
    ungroup() %>%
    group_by(term) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>% # p value correction
    filter(term %in% gradient_vec1)
  
  significant_unpropotional_025 <- clm_df_unpropotional %>% 
    dplyr::select(name,term,p.value) %>%
    group_by(term) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    group_by(name) %>%
    summarise(p.value = max(p.value),q.value = max(q.value)) %>%
    filter(p.value < 0.05, q.value < 0.25) %>%
    pluck("name")
  
  
  re_df <- data.frame(name = signif_assoc_cags)
  re_df_01 <- data.frame(name = signif_assoc_cags_01)
  re_df_025 <- data.frame(name = signif_assoc_cags_025)
  returnlist <- list(re_df=re_df,re_df_01 = re_df_01,re_df_025=re_df_025,clm_df=clm_df,
                     prop_odds_passed_cags = prop_odds_passed_cags,
                     clm_df_unpropotional = clm_df_unpropotional,
                     significant_unpropotional_025 = significant_unpropotional_025)
  return(returnlist)
}
Plot_heatmap <- function(df,vec,ord_df,name){
  # df <- Species_df
  # vec <- Species_list$re_df_025$name
  # ord_df <- Species_list$clm_df
  # name <- "Species"
  estimate_df <- ord_df %>%
    filter(name %in% vec) %>%
    mutate(sign = if_else(estimate > 0,"increase","decrease"))
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  Plot_df <- df %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble() %>%
    left_join(gradient) %>%
    relocate(c("Site","Season","Gradient"),.after = ID) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient, levels = c("High","Medium","Low"))) %>%
    arrange(Site) %>%
    dplyr::select(ID,Site,Season,Gradient,vec)
  Plot_mat <- Plot_df %>%
    dplyr::select(ID,vec) %>%
    column_to_rownames("ID") %>%
    as.matrix() %>%
    t()
  column_df <- data.frame(name = Plot_df$ID) %>%
    separate(name,c("Site","Season","tmp")) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient, levels = c("High","Medium","Low"))) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-tmp) %>%
    `rownames<-`(Plot_df$ID) %>%
    as.vector()
  bar <- columnAnnotation(df = column_df,col = list(Gradient = envPalette,
                                              Site = cbPalette,
                                              Season = seasonPalette),
                       #show_annotation_name = T,
                       show_legend = c(Gradient = FALSE),
                       annotation_legend_param = list(
                         Gradient = list(nrow = 1),
                         Site = list(nrow = 1),
                         Season = list(nrow = 1)))
  row_df <- estimate_df %>%
    dplyr::select(name,sign) %>%
    mutate(sign = factor(sign,levels = c("increase","decrease"))) %>%
    #arrange(sign) %>%
    column_to_rownames("name") %>%
    as.vector()
  pollutionPalette <- c("increase" ="#D55E00","decrease" = "#0072B2")
  foo <- rowAnnotation(df = row_df,col = list(sign = pollutionPalette),
                       show_annotation_name = FALSE,
                       show_legend = F,
                       annotation_legend_param = list(title = "",
                                                      nrow = 1))
  #library(circlize)
  #col_fun = colorRamp2(c(-4, 0, 4), c("#01befe", "white", "#ff7d00"))
  colors = structure(c("black","#ff7d00"), names = c("0", "1"))
  Ht = ComplexHeatmap::Heatmap(Plot_mat,cluster_columns = F,
                               #rect_gp = gpar(col = "white", lwd = 2),
                               top_annotation = bar,
                               right_annotation = foo,
                               column_split = column_df$Gradient,
                               cluster_column_slices = TRUE,
                               row_split = row_df$sign,
                               cluster_row_slices = TRUE,
                               col = colors,
                               row_names_max_width = max_text_width(rownames(Plot_mat),gp = gpar(fontsize = 12)),
                               column_names_max_height = max_text_width(colnames(Plot_mat),gp = gpar(fontsize = 12)),
                               heatmap_legend_param = list(title = name,nrow = 1,labels = c("absence","presence")))
  P <- draw(Ht,merge_legend = TRUE,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
}
## Results ---------------------------------------------------------------------
vec <- c("Low","Medium","High")
Species_list <- Ordinal_reg(Animal_species,vec)
return_list <- list(
  Species = Species_list
)
#load("~/Desktop/marine_sediment/results/Shelby_data/ordinal_reg_list.RData")
saveRDS(return_list, file="~/Desktop/marine_sediment/results_new/Shelby_data/ordinal_reg_list.rds")
# Phylum_P <- Plot_heatmap(Phylum_df,Phylum_list$re_df_025$name,Phylum_list$clm_df,"Phylum")
# Class_P <- Plot_heatmap(Class_df,Class_list$re_df_025$name,Class_list$clm_df,"Class")
# Order_P <- Plot_heatmap(Order_df,Order_list$re_df_025$name,Order_list$clm_df,"Order")
# Family_P <- Plot_heatmap(Family_df,Family_list$re_df_025$name,Family_list$clm_df,"Family")
# Genus_P <- Plot_heatmap(Genus_df,Genus_list$re_df_025$name,Genus_list$clm_df,"Genus")
Species_P <- Plot_heatmap(Animal_species,Species_list$re_df_01$name,Species_list$clm_df,"Species")
# 
# pdf("~/Desktop/marine_sediment/results/Shelby_data/Ordinal_reg_phylum.pdf",
#     width=13, height= 3)
# Phylum_P
# dev.off()
# pdf("~/Desktop/marine_sediment/results/Shelby_data/Ordinal_reg_class.pdf",
#     width=13, height= 4)
# Class_P
# dev.off()
# pdf("~/Desktop/marine_sediment/results/Shelby_data/Ordinal_reg_order.pdf",
#     width=13, height= 6)
# Order_P
# dev.off()
# pdf("~/Desktop/marine_sediment/results/Shelby_data/Ordinal_reg_family.pdf",
#     width=13, height= 11)
# Family_P
# dev.off()
# pdf("~/Desktop/marine_sediment/results/Shelby_data/Ordinal_reg_genus.pdf",
#     width=13, height= 12)
# Genus_P
# dev.off()
pdf("~/Desktop/marine_sediment/results_new/Shelby_data/Ordinal_reg_species.pdf",
    width=13, height= 12)
Species_P <- Plot_heatmap(Animal_species,Species_list$re_df_01$name,Species_list$clm_df,"Species")
dev.off()
