# Library  -----------------------------------------------------------------
#library(networkD3)
library(tidyverse)
library(dplyr)
library(plotly)
# ## Significant features --------------------------------------------------------
# data_path <- "~/Projects/Marine_Sediment/mediation/mediation_new/"
# 
# ancom2_gradient_list <- readRDS(paste0(data_path,"ancom2_gradient_list.rds"))
# 
# Microbiome_feature_tmp <- ancom2_gradient_list$Species$re_1$df %>%
#   filter(!W %in% "Inf") %>%
#   filter(detected_0.7 %in% "TRUE") %>% 
#   .$taxa_id
# Microbiome_feature <- ancom2_gradient_list$Species$re_1$res$eff_size %>% 
#   as.data.frame() %>%
#   filter(. < 0) %>%
#   rownames_to_column() %>%
#   filter(rowname %in% Microbiome_feature_tmp) %>%
#   .$rowname
# Microbiome_feature_pos_with_high_pollute <- setdiff(Microbiome_feature_tmp,Microbiome_feature)
# 
# Small_animal_org_reg_list <- readRDS(paste0(data_path,"ordinal_reg_list.rds"))
# Animal_TAXA <- read_csv(paste0(data_path,"TAXA_Linlin.csv"))
# colnames(Animal_TAXA)[1] <- "ID"
# Animal_TAXA <- Animal_TAXA %>%
#   dplyr::select(.,-ID,-Species) %>%
#   distinct()
# Small_animal_feature_genus <- Small_animal_org_reg_list$Genus$clm_df %>%
#   filter(estimate > 0) %>%
#   filter(name %in% Small_animal_org_reg_list$Genus$re_df_01$name) %>%
#   dplyr::select(name) %>%
#   left_join(Animal_TAXA,by = c("name" = "Genus"))
# 
# Small_animal_feature_tmp <- Small_animal_org_reg_list$Species$clm_df %>%
#   filter(estimate < 0) %>%
#   filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
# Small_animal_feature <- Small_animal_feature_tmp$name
# Small_animal_feature_tmp_pos <- Small_animal_org_reg_list$Species$clm_df %>%
#   filter(estimate > 0) %>%
#   filter(name %in% Small_animal_org_reg_list$Species$re_df_01$name)
# Small_animal_feature_pos <- Small_animal_feature_tmp_pos$name
# COG_list <- readRDS(paste0(data_path,"COG_ordinal_reg_list_1.rds"))
# COG_pathway_feature <- COG_list$COG_pathway_list$clm_df %>%
#   filter(estimate < 0) %>%
#   filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
# COG_pathway_feature_pos_with_high_pollute <- COG_list$COG_pathway_list$clm_df %>%
#   filter(estimate > 0) %>%
#   filter(name %in% COG_list$COG_pathway_list$re_df_01$name)
## Plot sankey -----------------------------------------------------------------
# re_df <- readRDS("~/Projects/Marine_Sediment/mediation/mediation_new/pollution_COG_benthos_neg.rds")
# re_df_tidy <- re_df %>%
#   mutate(BH_Beta_step2_iv = p.adjust(Beta_step2_iv_p)) %>%
#   mutate(BH_Beta_step3_med = p.adjust(Beta_step3_med_p)) %>%
#   #filter(BH_Beta_step2_iv < 0.25,BH_Beta_step3_med < 0.25) %>%
#   filter(!ACME_P_value %in% "") %>%
#   mutate(BH_P = p.adjust(ACME_P_value)) %>%
#   filter(ACME_P_value < 0.05,BH_P < 0.05) %>%
#   filter(Beta_step1_iv < 0,Beta_step3_iv < 0)

re_df_tidy <- read_tsv("~/Projects/Marine_Sediment/mediation/mediation_new/mediation_COGs.tsv")

tmp_df1 <- re_df_tidy %>%
  dplyr::select(iv,med) %>%
  mutate(name = iv) %>%
  `colnames<-`(c("source","target","name"))
tmp_df2 <- re_df_tidy %>%
  dplyr::select(med,dv,iv) %>%
  `colnames<-`(c("source","target","name"))

links <- rbind(tmp_df1,tmp_df2) #%>%
 # distinct()
Node <- unique(c(links$source,links$target))
IDsource <- match(links$source, Node)-1 
IDtarget <- match(links$target, Node)-1

Pollution_group <- read_tsv("~/Projects/Marine_Sediment/mediation/mediation/pollution_class.tsv")
Pollution_group$name <- str_remove_all(Pollution_group$name," \\(μg/kg\\)") %>%
  str_remove_all(.," \\(mg/kg\\)") %>%
  str_remove_all(.," \\(mg/L\\)") %>%
  str_remove_all(.," \\(μg/L\\)")
# Pollution_group <- Pollution_group %>%
#   filter(name %in% colnames(r)) %>%
#   mutate(name = factor(name,levels = colnames(r))) %>%
#   arrange(name) %>%
#   column_to_rownames("name")
Pollution_group_palette <- data.frame(
  group = c("Heavy metal","Nitrogen","PAH","Organic pollution"),
  color = c("#555555","#F3C57B","#9A8A76","#DB735C")
)
links_color <- links %>%
  left_join(Pollution_group) %>%
  left_join(Pollution_group_palette)

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
links_color1 <- links_color %>%
  group_by(color) %>%
  mutate(color1 = t_col(color))
 Color_df <- data.frame(name = Node) %>%
  left_join(Pollution_group) %>%
  left_join(Pollution_group_palette) %>%
  mutate(color = if_else(name %in% COG_pathway_feature$name,"#0072B2",color)) %>%
  mutate(color = if_else(name %in% COG_pathway_feature_pos_with_high_pollute$name,"#D55E00",color)) %>%
  mutate(color = if_else(name %in% Small_animal_feature,"#0072B2",color)) %>%
  mutate(color = if_else(name %in% Small_animal_feature_pos,"#D55E00",color)) %>%
  mutate(group = if_else(name %in% c(COG_pathway_feature$name,COG_pathway_feature_pos_with_high_pollute$name),"microbiome",group)) %>%
  mutate(group = if_else(name %in% c(Small_animal_feature,Small_animal_feature_pos),"Benthos",group))
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = Color_df$name,
    color = Color_df$color,
    #groups = Color_df$group,
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = IDsource,
    target = IDtarget,
    color = links_color1$color1,
    opacity = 0.5,
    value =  rep(1,150)
  )
)
fig <- fig %>% layout(
  #title = "Basic Sankey Diagram",
  font = list(
    size = 14
  )
)

fig
