# Title: Differential Abundance Plot
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Load library ----------------------------------------------------------------
library(tidyverse)
library(broom)
library(ComplexHeatmap)
library(ggsci)
## Load data -------------------------------------------------------------------
data_pollution <- readRDS("~/Desktop/marine_sediment/results/Benthos/Animal_List.rds")
data_season <- readRDS("~/Desktop/marine_sediment/results/Benthos/Animal_seasonList.rds")
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
  mutate(OTU = if_else(is.na(OTU) & !is.na(Phylum),paste0(Phylum," OTU",ID),OTU)) %>%
  dplyr::select(.,-ID) %>%
  distinct()
### Setting color parameters ---------------------------------------------------
envPalette <- c("High" ="#D55E00","Low" = "#0072B2")
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
seasonPalette <- c("Summer" = "#0073C2FF","Winter" = "#EFC000FF")
## Function --------------------------------------------------------------------
Plot_heatmap <- function(df_lst,taxa_file){
  # df_lst <- data_pollution
  # taxa_file <- Animal_TAXA

  estimate_df <- df_lst$Animal_otu_relist$sign_df %>%
    mutate(sign = if_else(coef > 0,"increase","decrease")) %>%
    dplyr::select(feature,sign) %>%
    left_join(df_lst$Animal_otu_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  
  Plot_df <- df_lst$Animal_otu_list$taxon %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    pivot_longer(-1) %>%
    filter(name %in% estimate_df$feature) %>%
    left_join(estimate_df %>% dplyr::rename(name = feature)) %>%
    dplyr::select(.,-c(name,sign)) %>%
    dplyr::rename(name = Taxa) %>%
    pivot_wider() %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble() %>%
    left_join(gradient) %>%
    relocate(c("Site","Season","Gradient"),.after = ID) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient, levels = c("High","Low"))) %>%
    arrange(Site)
  
  Plot_mat <- Plot_df %>%
    dplyr::select(.,-c(Site,Season,Gradient)) %>%
    column_to_rownames("ID") %>%
    as.matrix() %>%
    t()
  
  all(colnames(Plot_mat) == (Plot_df$ID))
  
  column_df <- data.frame(name = colnames(Plot_mat)) %>%
    separate(name,c("Site","Season","tmp")) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    dplyr::rename(Pollution = Gradient) %>%
    mutate(Pollution = factor(Pollution, levels = c("High","Low"))) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-tmp) %>%
    `rownames<-`(colnames(Plot_mat)) %>%
    as.vector()
  
  all(colnames(Plot_mat) == rownames(column_df))
  
  bar <- columnAnnotation(df = column_df,col = list(Pollution = envPalette,
                                              Site = cbPalette,
                                              Season = seasonPalette),
                       #show_annotation_name = T,
                       #show_legend = c(Pollution = FALSE),
                       annotation_legend_param = list(
                         Pollution = list(nrow = 1,
                                          title_gp = gpar(fontsize = 16, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 16, fontface = "bold")
                                          ),
                         Site = list(nrow = 1,
                                     title_gp = gpar(fontsize = 16, fontface = "bold"),
                                     labels_gp = gpar(fontsize = 16, fontface = "bold")
                                     ),
                         Season = list(nrow = 1,
                                       title_gp = gpar(fontsize = 16, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 16, fontface = "bold")
                                       )#,
                         #title_gp = gpar(fontsize = 16, fontface = "bold"),  # Increase font size of the legend title
                         #labels_gp = gpar(fontsize = 16, fontface = "bold")
                         ),
                       annotation_name_gp= gpar(fontsize = 16, fontface = "bold")
                       )
  
  row_df <- data.frame(Taxa = rownames(Plot_mat)) %>%
    left_join(estimate_df) %>%
    dplyr::rename(OTU = Taxa) %>%
    left_join(taxa_file) %>%
    dplyr::select(OTU,Phylum,sign) %>%
    mutate(sign = factor(sign,levels = c("increase","decrease"))) %>%
    #arrange(sign) %>%
    column_to_rownames("OTU") %>%
    as.vector()
  
  Phylum_names <- unique(row_df$Phylum)
  set.seed(123)
  P10 = pal_simpsons("springfield")(11)
  #swatch(P28)
  names(P10) <- Phylum_names
  saveRDS(P10,"~/Desktop/marine_sediment/results/Benthos/Animal_phylum_color.rds")
  
  pollutionPalette <- c("increase" ="#D55E00","decrease" = "#0072B2")
  foo <- rowAnnotation(df = row_df,col = list(sign = pollutionPalette,
                                              Phylum = P10),
                       #show_annotation_name = FALSE,
                       show_legend = c(sign = FALSE),
                       annotation_legend_param = list(
                         sign = list(nrow = 1,
                                     title_gp = gpar(fontsize = 16, fontface = "bold"),
                                     labels_gp = gpar(fontsize = 16, fontface = "bold")
                                     ),
                         Phylum = list(nrow = 1,
                                       title_gp = gpar(fontsize = 16, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 16, fontface = "bold")
                                       )
                       ),
                       annotation_name_gp= gpar(fontsize = 16, fontface = "bold")
                       )
  #library(circlize)
  #col_fun = colorRamp2(c(-4, 0, 4), c("#01befe", "white", "#ff7d00"))
  colors = structure(c("black","#ff7d00"), names = c("0", "1"))
  Ht = ComplexHeatmap::Heatmap(Plot_mat,cluster_columns = F,
                               #rect_gp = gpar(col = "white", lwd = 2),
                               top_annotation = bar,
                               right_annotation = foo,
                               column_split = column_df$Pollution,
                               cluster_column_slices = TRUE,
                               row_split = row_df$sign,
                               cluster_row_slices = TRUE,
                               col = colors,
                               column_names_gp = gpar(fontsize = 16, fontface = "bold"),
                               row_names_gp = gpar(fontsize = 16, fontface = "bold"),
                               column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                               row_title_gp = gpar(fontsize = 16, fontface = "bold"),
                               row_names_max_width = max_text_width(rownames(Plot_mat),gp = gpar(fontsize = 16, fontface = "bold")),
                               column_names_max_height = max_text_width(colnames(Plot_mat),gp = gpar(fontsize = 16, fontface = "bold")),
                               heatmap_legend_param = list(title = "Presence/Absence",
                                                           nrow = 1,labels = c("absence","presence"),
                                                           title_gp = gpar(fontsize = 16, fontface = "bold"),  # Adjust the legend title font size
                                                           labels_gp = gpar(fontsize = 16, fontface = "bold")
                                                           )
                               )
  P <- draw(Ht,merge_legend = TRUE,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")
  return(Ht)
}
Plot_heatmap_season <- function(df_lst,taxa_file){
  # df_lst <- data_season
  # taxa_file <- Animal_TAXA
  
  estimate_df <- df_lst$Animal_otu_relist$sign_df %>%
    mutate(sign = if_else(coef > 0,"increase","decrease")) %>%
    dplyr::select(feature,sign) %>%
    left_join(df_lst$Animal_otu_list$Feature_ID %>% dplyr::rename(feature = ID,Taxa = name))
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  
  Plot_df <- df_lst$Animal_otu_list$taxon %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    pivot_longer(-1) %>%
    filter(name %in% estimate_df$feature) %>%
    left_join(estimate_df %>% dplyr::rename(name = feature)) %>%
    dplyr::select(.,-c(name,sign)) %>%
    dplyr::rename(name = Taxa) %>%
    pivot_wider() %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble() %>%
    left_join(gradient) %>%
    relocate(c("Site","Season","Gradient"),.after = ID) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient, levels = c("High","Low"))) %>%
    arrange(Site)
  
  Plot_mat <- Plot_df %>%
    dplyr::select(.,-c(Site,Season,Gradient)) %>%
    column_to_rownames("ID") %>%
    as.matrix() %>%
    t()
  
  all(colnames(Plot_mat) == (Plot_df$ID))
  
  column_df <- data.frame(name = colnames(Plot_mat)) %>%
    separate(name,c("Site","Season","tmp")) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    dplyr::rename(Pollution = Gradient) %>%
    mutate(Pollution = factor(Pollution, levels = c("High","Low"))) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-tmp) %>%
    `rownames<-`(colnames(Plot_mat)) %>%
    as.vector()
  
  all(colnames(Plot_mat) == rownames(column_df))
  
  bar <- columnAnnotation(df = column_df,col = list(Pollution = envPalette,
                                                    Site = cbPalette,
                                                    Season = seasonPalette),
                          #show_annotation_name = T,
                          show_legend = c(Pollution = FALSE),
                          annotation_legend_param = list(
                            Pollution = list(nrow = 1),
                            Site = list(nrow = 1),
                            Season = list(nrow = 1)))
  
  row_df <- data.frame(Taxa = rownames(Plot_mat)) %>%
    left_join(estimate_df) %>%
    dplyr::rename(OTU = Taxa) %>%
    left_join(taxa_file) %>%
    dplyr::select(OTU,Phylum,sign) %>%
    mutate(sign = factor(sign,levels = c("increase","decrease"))) %>%
    #arrange(sign) %>%
    column_to_rownames("OTU") %>%
    as.vector()
  
  # Phylum_names <- unique(row_df$Phylum)
  # set.seed(123)
  # P10 = pal_simpsons("springfield")(11)
  # #swatch(P28)
  # names(P10) <- Phylum_names
  P10 <- readRDS("~/Desktop/marine_sediment/results/Benthos/Animal_phylum_color.rds")
  
  pollutionPalette <- c("increase" ="#EFC000FF","decrease" = "#0073C2FF")
  foo <- rowAnnotation(df = row_df,col = list(sign = pollutionPalette,
                                              Phylum = P10),
                       #show_annotation_name = FALSE,
                       show_legend = c(sign = FALSE),
                       annotation_legend_param = list(
                         sign = list(nrow = 1),
                         Phylum = list(nrow = 1)
                       ))
  #library(circlize)
  #col_fun = colorRamp2(c(-4, 0, 4), c("#01befe", "white", "#ff7d00"))
  colors = structure(c("black","#ff7d00"), names = c("0", "1"))
  Ht = ComplexHeatmap::Heatmap(Plot_mat,cluster_columns = F,
                               #rect_gp = gpar(col = "white", lwd = 2),
                               top_annotation = bar,
                               right_annotation = foo,
                               column_split = column_df$Season,
                               cluster_column_slices = TRUE,
                               row_split = row_df$sign,
                               cluster_row_slices = TRUE,
                               col = colors,
                               row_names_max_width = max_text_width(rownames(Plot_mat),gp = gpar(fontsize = 12)),
                               column_names_max_height = max_text_width(colnames(Plot_mat),gp = gpar(fontsize = 12)),
                               heatmap_legend_param = list(title = "Presence/Absence",nrow = 1,labels = c("absence","presence")))
  P <- draw(Ht,merge_legend = TRUE,
            heatmap_legend_side = "bottom",
            annotation_legend_side = "bottom")
  return(Ht)
}

## Results ---------------------------------------------------------------------
Species_P <- Plot_heatmap(data_pollution,Animal_TAXA)

P <- draw(Species_P,merge_legend = TRUE,
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom")

pdf("~/Desktop/marine_sediment/results/Benthos/DA_pollution_plot.pdf",
    width=22, height= 25)
draw(Species_P,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

## Result season ---------------------------------------------------------------
Species_P <- Plot_heatmap_season(data_season,Animal_TAXA)

P <- draw(Species_P,merge_legend = TRUE,
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom")

pdf("~/Desktop/marine_sediment/results/Benthos/DA_season_plot.pdf",
    width=13, height= 16)
draw(Species_P,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()
