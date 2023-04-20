## Load library ----------------------------------------------------------------
library(ordinal)
library(tidyverse)
library(broom)
library(ComplexHeatmap)
library(janitor)
library(ggprism)
set.seed(123)
data_path <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/"
COG_pathway <- read_tsv(paste0(data_path,"COG_pathway_stat.flt.hellinger.tsv"))
COG <- read_tsv(paste0(data_path,"COG_cog.flt.hellinger.tsv"))
COG_cate <- read_tsv(paste0(data_path,"COG_cate.tsv")) %>%
  column_to_rownames("Category") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column()
### Setting color parameters ---------------------------------------------------
cbPalette1 <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
envPalette1 <- c("#0072B2","#E69F00","#D55E00")
envPalette <- c("High" ="#D55E00","Medium" = "#E69F00","Low" = "#0072B2")
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
seasonPalette <- c("Summer" = "#0073C2FF","Winter" = "#EFC000FF")
gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                   Gradient = c("High","High","High","Medium","Medium","Low","Low"))
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
  #df <- COG
  #vec <- c("Low","Medium","High")
  ## Tidy data -----------------------------------------------------------------
  taxon <- df %>%
    replace(is.na(.),0) %>%
    column_to_rownames()
  taxon[1:length(taxon)] <- apply(taxon[1:length(taxon)],2,as.numeric)
  site_name <- sapply(strsplit(rownames(taxon),'_'), function(v) return(v[1]))
  season_name <- sapply(strsplit(rownames(taxon),'_'), function(v) return(v[2])) %>%
    gsub("S","Summer",.) %>% gsub("W","Winter",.)
  ## 50% prevalence filter in at leat one site
  Prev_df <- cbind(taxon,Site =site_name) %>%
    rownames_to_column() %>%
    pivot_longer(2:(length(.)-1),names_to = "microbiome",values_to = "relativeAbundance") %>%
    group_by(microbiome,Site) %>%
    summarise_at("relativeAbundance", function(x) length(which(x!=0))) %>%
    ungroup() %>%
    group_by(microbiome) %>%
    summarise_at("relativeAbundance", function(x) length(which(x>=3))) %>%
    filter(relativeAbundance > 0)
  taxon <- taxon %>%
    rownames_to_column() %>%
    pivot_longer(2:length(.),names_to = "microbiome",values_to = "relativeAbundance") %>%
    filter(microbiome %in% Prev_df$microbiome) %>%
    pivot_wider(names_from = microbiome,values_from = relativeAbundance) %>%
    column_to_rownames() #%>%
    # decostand("hellinger")
  
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  df <- taxon %>%
    rownames_to_column("ID") %>%
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
  
  
  re_df <- data.frame(name = signif_assoc_cags)
  re_df_01 <- data.frame(name = signif_assoc_cags_01)
  re_df_025 <- data.frame(name = signif_assoc_cags_025)
  returnlist <- list(re_df=re_df,re_df_025=re_df_025,re_df_01 = re_df_01,clm_df=clm_df,
                     prop_odds_passed_cags = prop_odds_passed_cags)
  return(returnlist)
}
Plot_gradient <- function(taxon,feature_df){
  # taxon <- COG_cate
  # feature_df <- COG_cate_list$re_df_01$name
  feature <- feature_df
  RelativeAbundance_df <- taxon %>%
    rename(SampleID = rowname) %>%
    mutate(tmp = SampleID) %>%
    separate(tmp,c("Site", "Season","tmp")) %>%
    dplyr::select(.,-tmp) %>%
    relocate(c("Site", "Season"),.after = SampleID) %>%
    as_tibble() %>%
    mutate_all(type.convert) %>%
    pivot_longer(4:(length(.)),names_to = "microbiome",values_to = "relativeAbundance") %>%
    filter(microbiome %in% feature) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient,levels = c("Low","Medium","High")))
  p <- ggplot(RelativeAbundance_df) + 
    geom_boxplot(aes(x = Gradient, y = relativeAbundance,fill = Gradient),alpha = 0.2) +
    facet_wrap(~ microbiome,scales = "free",ncol = 6) +
    geom_point(aes(x = Gradient, y = relativeAbundance,color = Site,shape = Season),size=1.5) +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = envPalette)+
    scale_shape_manual(values = c(17,16)) +
    xlab("") + ylab("")+ theme_classic()+ theme(legend.title = element_blank()) +
    theme(strip.background = element_blank()) +
    theme_prism(base_size = 16) +
    theme(legend.position="none")
  return(p)
}
## Results ---------------------------------------------------------------------
vec <- c("Low","Medium","High")
COG_pathway_list <- Ordinal_reg(COG_pathway,vec)
COG_list <- Ordinal_reg(COG,vec)
COG_cate_list <- Ordinal_reg(COG_cate,vec)
  
return_list <- list(
  COG_pathway_list = COG_pathway_list,
  COG_list = COG_list,
  COG_cate_list = COG_cate_list
)
saveRDS(return_list, file="~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_ordinal_reg_list_1.rds")
return_list <- readRDS("~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_ordinal_reg_list_1.rds")
COG_pathway_P <- Plot_gradient(COG_pathway,return_list$COG_pathway_list$re_df_01$name)
COG_cate_p <- Plot_gradient(COG_cate,return_list$COG_cate_list$re_df_01$name)
COG_p <- Plot_gradient(COG,c("COG0005","COG3096"))
#COG0005 increase from low to high
pdf("~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiome_v3/Figures_v3/functional_analysesOrdinal_reg_cog_cate_new.pdf",
    width=18, height= 5)
COG_cate_p
dev.off()

pdf("~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiome_v5/Marine_microbiome_v5/Figures/functional_cog_pathway.pdf",
    width=25, height= 18)
COG_pathway_P
dev.off()

write_tsv(return_list$COG_pathway_list$re_df_01,"~/Desktop/marine_sediment/results_new/functional_analyses/COG_pathway_1.tsv")
pdf("~/Desktop/marine_sediment/results_new/functional_analysesOrdinal_reg_cog_cate_new.pdf",
    width=10, height= 5)
COG_cate_p
dev.off()
COG_df <- return_list$COG_list$clm_df %>%
  filter(name %in% return_list$COG_list$re_df_01$name) %>%
  mutate(sign = if_else(estimate < 0,"neg","pos"))
cate_df <- read_tsv("~/Desktop/marine_sediment/results_new/functional_analyses/fun-20.tab")
cate_tab <- COG_cate_list$re_df_01 %>%
  rename(ID = "name") %>%
  left_join(cate_df)
write_tsv("~/Desktop/marine_sediment/results_new/functional_analyses/sign_cate.tsv")
COG_df <- return_list$COG_list$clm_df %>%
  filter(name %in% return_list$COG_list$re_df_01$name)
write_tsv(COG_df,"~/Desktop/marine_sediment/manuscript/Figure_v5/TabS8_COG_items.tsv")
