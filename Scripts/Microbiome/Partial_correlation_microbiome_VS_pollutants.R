## Load library ----------------------------------------------------------------
library(tidyverse)
library(Hmisc) ## calculate microbiom correlation coefficient
library(psych)
library(spearmanCI)
library(vegan)
library('RColorBrewer')
library(ComplexHeatmap)
library(ppcor)
data_path <- "~/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
## pollution data: pollution ---------------------------------------------------
## load sediment data -------------------------------------------------------------------
file1 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ms5.csv")
file2 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ns6.csv")
file3 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ps6.csv")
file4 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ss1-ss3-ss5.csv")
file5 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ts2.csv")
file6 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5yrs-ms8.csv")
file_sediment <- rbind(file1,file2,file3,file4,file5,file6)
## pollution data: pollution ---------------------------------------------------
pollution_sediment <- read_csv2("~/Desktop/marine_sediment/data/pollution.csv") %>%
  dplyr::select(STATION,StationID) %>%
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
  dplyr::rename(Site = StationID, Season = season) %>%
  type_convert() %>%
  dplyr::select(where(~ !is.numeric(.) || var(., na.rm = TRUE) != 0)) %>%
  pivot_longer(3:length(.)) %>%
  mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
  mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
  dplyr::rename(meta_name = name, meta_value = value)
## Load data -------------------------------------------------------------------
gradient_list <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/Microbiome_DA_pollution_List.rds")
taxa_tab <- read_tsv("~/Desktop/marine_sediment/MEGAN6/STAMP_out/taxa.tsv") %>%
  dplyr::select(.,-IDs) %>%
  distinct()
## Partial Spearman correlation ------------------------------------------------
Partial_Spearman_adjust_season <- function(data_list,feature_list,meta_data)
{
  # data_list <- gradient_list$species_list
  # feature_list <- gradient_list$species_relist
  # meta_data <- Score_sediment
  
  pcor_test <- function(df,x,y,p){
    a <- df[[x]]
    b <- df[[y]]
    c <- df[[p]] %>% as.factor() %>% as.integer()
    r <- as_tibble(pcor.test(a,b,c,method = "spearman"))
    return(r)
  }
  
  microbe <- data_list$taxon %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    pivot_longer(-1, names_to = "ID") %>%
    filter(ID %in% feature_list$sign_df$feature) %>%
    left_join(data_list$Feature_ID) %>%
    left_join(data_list$meta_data %>% rownames_to_column()) %>%
    group_by(Gradient,Site,Season,name) %>%
    summarise_at("value", mean, na.rm = TRUE) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient,levels = c("Low","High"))) %>%
    left_join(meta_data) %>%
    group_by(name,meta_name) %>%
    group_modify(~ pcor_test(.x,"value","meta_value","Season")) %>%
    bind_rows() %>%
    ungroup() %>%
    group_by(meta_name) %>%
    mutate(q.value = p.adjust(p.value, method = "BH")) #%>% # p value correction
  
  r <- microbe %>%
    dplyr::select(name,meta_name,estimate) %>%
    pivot_wider(names_from = meta_name,values_from = estimate) %>%
    column_to_rownames(var = "name") %>%
    as.matrix()
  p <- microbe %>%
    dplyr::select(name,meta_name,p.value) %>%
    pivot_wider(names_from = meta_name,values_from = p.value) %>%
    column_to_rownames(var = "name") %>%
    as.matrix()
  p.BH <- microbe %>%
    dplyr::select(name,meta_name,q.value) %>%
    pivot_wider(names_from = meta_name,values_from = q.value) %>%
    column_to_rownames(var = "name") %>%
    as.matrix()
  p_1 <- p
  p_1[p >= 0.05] <- 1
  p_1[p < 0.05] <- 0
  p.BH <- p.BH + p_1
  ## Heatmap ---------------------------------------------------------------------
  re_p <- Heatmap(r,name = "Spearman r",rect_gp = gpar(col = "white", lwd = 2),
          row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 12)),
          column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 12)),
          layer_fun = function(j, i, x, y, width, height, fill) {
            # here i is the row number from 1-6 and from 7-13
            # j is the the column number from 1-57 ( two list : the top part and the down part)
            # x coordinate of middle point of the cell which is measured in the viewport of the heatmap body.
            # y coordinate of middle point of the cell which is measured in the viewport of the heatmap body
            ind_mat = restore_matrix(j, i, x, y)  
            
            
            v = pindex(p, i, j)
            l = v < 0.05
            ind_mat_vec = as.matrix(ind_mat)
            ind = ind_mat_vec[l]
            grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
            
            
            v2 = pindex(p.BH, i, j)
            l2 = v2 < 0.2
            ind_mat_vec2 = as.matrix(ind_mat)
            ind2 = ind_mat_vec2[l2]
            grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))

          }
          )
    return_list <- list(df = microbe,
                        meta_data = Score_sediment,
                        p = re_p)
    return(return_list)
}
## Results ---------------------------------------------------------------------
Phylum_list <- Partial_Spearman_adjust_season(gradient_list$phylum_list,gradient_list$phylum_relist,Score_sediment)
Genus_list <- Partial_Spearman_adjust_season(gradient_list$genus_list,gradient_list$genus_relist,Score_sediment)
Species_list <- Partial_Spearman_adjust_season(gradient_list$species_list,gradient_list$species_relist,Score_sediment)
## Figure in manuscript --------------------------------------------------------
sign_spe <- gradient_list$species_relist$sign_df %>%
  dplyr::rename(ID = feature,tmp = name) %>%
  left_join(gradient_list$species_list$Feature_ID) %>%
  dplyr::rename(Species = name) %>%
  left_join(taxa_tab) %>%
  filter(Phylum %in% c("p__Proteobacteria","p__Planctomycetes",
                       "p__Nitrospinae","p__Chloroflexi",
                       "p__Candidatus Tectomicrobia","p__Bacteroidetes"
  )
  ) %>%
  filter(!grepl("^g__\\(", Genus)) %>%
  filter(!grepl("sp\\.", Species)) %>%
  dplyr::rename(feature = ID)

sign_genus <- gradient_list$genus_relist$sign_df %>%
  dplyr::rename(ID = feature,tmp = name) %>%
  left_join(gradient_list$genus_list$Feature_ID) %>%
  dplyr::rename(Genus = name) %>%
  left_join(taxa_tab %>% dplyr::select(.,-Species) %>% distinct()) %>%
  filter(Genus %in% sign_spe$Genus) %>%
  dplyr::rename(feature = ID)

sign_phylum <- gradient_list$phylum_relist$sign_df %>%
  dplyr::rename(ID = feature,tmp = name) %>%
  left_join(gradient_list$phylum_list$Feature_ID) %>%
  filter(name %in% c("p__Proteobacteria","p__Planctomycetes",
                     "p__Nitrospinae","p__Chloroflexi",
                     "p__Candidatus Tectomicrobia","p__Bacteroidetes"
  )) %>%
  dplyr::rename(feature = ID)

Sign_Species_lst <- list(sign_df = sign_spe)
Sign_Genus_lst <- list(sign_df = sign_genus)
Sign_Phylum_lst <- list(sign_df = sign_phylum)
## results in manuscript -------------------------------------------------------
Phylum_manu_list <- Partial_Spearman_adjust_season(gradient_list$phylum_list,Sign_Phylum_lst,Score_sediment)
Genus_manu_list <- Partial_Spearman_adjust_season(gradient_list$genus_list,Sign_Genus_lst,Score_sediment)
Species_manu_list <- Partial_Spearman_adjust_season(gradient_list$species_list,Sign_Species_lst,Score_sediment)
## Heatmap prepare -------------------------------------------------------------
sign_tab <- rbind(
  data.frame(Phylum_manu_list$df,level="Phylum"),
  data.frame(Genus_manu_list$df,level="Genus"),
  data.frame(Species_manu_list$df,level="Species")
) %>%
  mutate(name = str_remove_all(name,"^.*?__"))
level_vec <- c("Phylum","Genus","Species")
cbPalette <- c("Phylum" = "#E69F00","Genus"= "#0072B2","Species" = "#D55E00")
r <- sign_tab %>%
  dplyr::select(name,meta_name,estimate) %>%
  pivot_wider(names_from = meta_name,values_from = estimate) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p <- sign_tab %>%
  dplyr::select(name,meta_name,p.value) %>%
  pivot_wider(names_from = meta_name,values_from = p.value) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p.BH <- sign_tab %>%
  dplyr::select(name,meta_name,q.value) %>%
  pivot_wider(names_from = meta_name,values_from = q.value) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p_1 <- p
p_1[p >= 0.05] <- 1
p_1[p < 0.05] <- 0
p.BH <- p.BH + p_1

## Heatmap ---------------------------------------------------------------------
Pollution_group <- read_tsv("~/Desktop/marine_sediment/renew_pollutionbin/sediment_pollute.tsv")
Pollution_group <- Pollution_group %>%
  mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
  mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
  dplyr::select(name,categories)
Pollution_group <- Pollution_group %>%
  filter(name %in% colnames(r)) %>%
  mutate(name = factor(name,levels = colnames(r))) %>%
  arrange(name) %>%
  column_to_rownames("name")
all_sign_plot <- sign_tab %>% ungroup() %>%  filter(meta_name == "Lead") %>%
  mutate(level = factor(level,levels = level_vec))
all(all_sign_plot$name == rownames(r))
rownames(all_sign_plot) <- all_sign_plot$name
Sign = data.frame(level = all_sign_plot$level)
#cbPalette <- c("Genus"= "#0072B2","Species" = "#D55E00")
bar <- columnAnnotation(df = Sign,col = list(level =cbPalette),
                        annotation_name_gp = gpar(fontsize = 16),
                        annotation_legend_param = list(nrow = 1,
                                                       labels_gp = gpar(fontsize = 16),
                                                       title_gp = gpar(fontsize = 16,fontface = "bold"))
)
pollutionPalette <- c("Heavy Metal" ="#555555",
                      "Metal" = "grey",
                      "Inorganic pollutants" = "#F3C57B",
                      "PAH" = "#9A8A76",
                      "Nutrient pollution"="#DB735C")
foo <- rowAnnotation(df = Pollution_group,col = list(categories = pollutionPalette),
                     annotation_name_gp = gpar(fontsize = 16),
                     annotation_legend_param = list(title = "Pollution group",
                                                    nrow = 1,
                                                    labels_gp = gpar(fontsize = 16),
                                                    title_gp = gpar(fontsize = 16,fontface = "bold")))

ht_list = Heatmap(t(r),name = "Partial Spearman r",rect_gp = gpar(col = "white", lwd = 2),
                  row_names_max_width = max_text_width(rownames(t(r)),gp = gpar(fontsize = 16)),
                  column_names_max_height = max_text_width(colnames(t(r)),gp = gpar(fontsize = 16)),
                  #cluster_columns = FALSE,
                  right_annotation = foo,
                  column_names_rot = 75,
                  column_names_gp = grid::gpar(fontsize = 16),
                  row_names_gp = grid::gpar(fontsize = 16),
                  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                  column_split = factor(Sign$level,levels = c("Phylum","Genus","Species")),
                  cluster_column_slices = F,
                  row_title = NULL,
                  #row_split = Sign$level,
                  #cluster_columns = F,
                  # label the circle of fdr 
                  #layer_fun = function(j, i, x, y, width, height, fill) {
                  #  ind_mat = restore_matrix(j, i, x, y)
                  #  v = pindex(p.BH, i, j) # qvalue is the matrix : qvalue table 
                  #  l = v < 0.2
                  #  ind_mat_vec = as.matrix(ind_mat)
                  #  ind = ind_mat_vec[l]
                  #  grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                  #},
                  top_annotation = bar,
                  heatmap_legend_param = list(direction = "horizontal",
                                              labels_gp = gpar(fontsize = 16),
                                              title_gp = gpar(fontsize = 16,fontface = "bold")),
                  #bottom_annotation = bar,
                  #cell_fun = function(j, i, x, y, width, height, fill) {
                  #  if(p[i, j] < 0.05)
                  #    grid.text("*", x, y, gp = gpar(fontsize = 10))
                  #},
                  layer_fun = function(j, i, x, y, width, height, fill) {
                    # here i is the row number from 1-6 and from 7-13
                    # j is the the column number from 1-57 ( two list : the top part and the down part)
                    # x coordinate of middle point of the cell which is measured in the viewport of the heatmap body.
                    # y coordinate of middle point of the cell which is measured in the viewport of the heatmap body
                    ind_mat = restore_matrix(j, i, x, y)  
                    
                    
                    v = pindex(t(p), i, j)
                    l = v < 0.05
                    ind_mat_vec = as.matrix(ind_mat)
                    ind = ind_mat_vec[l]
                    if (length(ind) != 0){
                      grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                    }
                    #grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                    
                    
                    v2 = pindex(t(p.BH), i, j)
                    l2 = v2 < 0.2
                    ind_mat_vec2 = as.matrix(ind_mat)
                    ind2 = ind_mat_vec2[l2]
                    if (length(ind2) != 0){
                      grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))
                    }
                    #grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))
                    
                  })
#dev.off()
pdf("~/Desktop/marine_sediment/results/MEGAN6/Correlation_microbiome_pollutants.pdf",
    height = 12,width = 23)
draw(ht_list,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()



## Heatmap prepare1 -------------------------------------------------------------
sign_tab <- rbind(
  data.frame(Phylum_list$df,level="Phylum"),
  data.frame(Genus_list$df,level="Genus"),
  data.frame(Species_list$df,level="Species")
) %>%
  mutate(name = str_remove_all(name,"^.*?__"))
level_vec <- c("Phylum","Genus","Species")
cbPalette <- c("Phylum" = "#E69F00","Genus"= "#0072B2","Species" = "#D55E00")
r <- sign_tab %>%
  dplyr::select(name,meta_name,estimate) %>%
  pivot_wider(names_from = meta_name,values_from = estimate) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p <- sign_tab %>%
  dplyr::select(name,meta_name,p.value) %>%
  pivot_wider(names_from = meta_name,values_from = p.value) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p.BH <- sign_tab %>%
  dplyr::select(name,meta_name,q.value) %>%
  pivot_wider(names_from = meta_name,values_from = q.value) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()
p_1 <- p
p_1[p >= 0.05] <- 1
p_1[p < 0.05] <- 0
p.BH <- p.BH + p_1

## Heatmap ---------------------------------------------------------------------
Pollution_group <- read_tsv("~/Desktop/marine_sediment/renew_pollutionbin/sediment_pollute.tsv")
Pollution_group <- Pollution_group %>%
  mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
  mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
  dplyr::select(name,categories)
Pollution_group <- Pollution_group %>%
  filter(name %in% colnames(r)) %>%
  mutate(name = factor(name,levels = colnames(r))) %>%
  arrange(name) %>%
  column_to_rownames("name")
all_sign_plot <- sign_tab %>% ungroup() %>%  filter(meta_name == "Lead") %>%
  mutate(level = factor(level,levels = level_vec))
all(all_sign_plot$name == rownames(r))
rownames(all_sign_plot) <- all_sign_plot$name
Sign = data.frame(level = all_sign_plot$level)
#cbPalette <- c("Genus"= "#0072B2","Species" = "#D55E00")
bar <- rowAnnotation(df = Sign,col = list(level =cbPalette),
                        annotation_name_gp = gpar(fontsize = 16),
                        annotation_legend_param = list(nrow = 1,
                                                       labels_gp = gpar(fontsize = 16),
                                                       title_gp = gpar(fontsize = 16,fontface = "bold"))
)
pollutionPalette <- c("Heavy Metal" ="#555555",
                      "Metal" = "grey",
                      "Inorganic pollutants" = "#F3C57B",
                      "PAH" = "#9A8A76",
                      "Nutrient pollution"="#DB735C")
foo <- columnAnnotation(df = Pollution_group,col = list(categories = pollutionPalette),
                     annotation_name_gp = gpar(fontsize = 16),
                     annotation_legend_param = list(title = "Pollution group",
                                                    nrow = 1,
                                                    labels_gp = gpar(fontsize = 16),
                                                    title_gp = gpar(fontsize = 16,fontface = "bold")))

ht_list = Heatmap(r,name = "Partial Spearman r",rect_gp = gpar(col = "white", lwd = 2),
                  row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 16)),
                  column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 16)),
                  #cluster_columns = FALSE,
                  right_annotation = bar,
                  column_names_rot = 75,
                  column_names_gp = grid::gpar(fontsize = 16),
                  row_names_gp = grid::gpar(fontsize = 16),
                  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                  row_split = factor(Sign$level,levels = c("Phylum","Genus","Species")),
                  cluster_row_slices = F,
                  row_title = NULL,
                  #row_split = Sign$level,
                  #cluster_columns = F,
                  # label the circle of fdr 
                  #layer_fun = function(j, i, x, y, width, height, fill) {
                  #  ind_mat = restore_matrix(j, i, x, y)
                  #  v = pindex(p.BH, i, j) # qvalue is the matrix : qvalue table 
                  #  l = v < 0.2
                  #  ind_mat_vec = as.matrix(ind_mat)
                  #  ind = ind_mat_vec[l]
                  #  grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                  #},
                  top_annotation = foo,
                  heatmap_legend_param = list(direction = "horizontal",
                                              labels_gp = gpar(fontsize = 16),
                                              title_gp = gpar(fontsize = 16,fontface = "bold")),
                  #bottom_annotation = bar,
                  #cell_fun = function(j, i, x, y, width, height, fill) {
                  #  if(p[i, j] < 0.05)
                  #    grid.text("*", x, y, gp = gpar(fontsize = 10))
                  #},
                  layer_fun = function(j, i, x, y, width, height, fill) {
                    # here i is the row number from 1-6 and from 7-13
                    # j is the the column number from 1-57 ( two list : the top part and the down part)
                    # x coordinate of middle point of the cell which is measured in the viewport of the heatmap body.
                    # y coordinate of middle point of the cell which is measured in the viewport of the heatmap body
                    ind_mat = restore_matrix(j, i, x, y)  
                    
                    
                    v = pindex(p, i, j)
                    l = v < 0.05
                    ind_mat_vec = as.matrix(ind_mat)
                    ind = ind_mat_vec[l]
                    if (length(ind) != 0){
                      grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                    }
                    #grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
                    
                    
                    v2 = pindex(p.BH, i, j)
                    l2 = v2 < 0.2
                    ind_mat_vec2 = as.matrix(ind_mat)
                    ind2 = ind_mat_vec2[l2]
                    if (length(ind2) != 0){
                      grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))
                    }
                    #grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))
                    
                  })
#dev.off()
pdf("~/Desktop/marine_sediment/results/MEGAN6/Correlation_microbiome_pollutants_supp.pdf",
    height = 52,width = 23)
draw(ht_list,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()


