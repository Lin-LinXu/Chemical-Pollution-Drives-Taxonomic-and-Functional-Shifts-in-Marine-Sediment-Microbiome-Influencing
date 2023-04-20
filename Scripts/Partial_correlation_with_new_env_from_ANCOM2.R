Partial_correlation_with_env_from_ANCOM2 <- function()
  {
## Load library ----------------------------------------------------------------
library(tidyverse)
library(Hmisc) ## calculate microbiom correlation coefficient
library(psych)
library(spearmanCI)
library(vegan)
library('RColorBrewer')
library(ComplexHeatmap)
library(ppcor)
data_path <- "~/Desktop/marine_sediment/data/"
## pollution data: pollution ---------------------------------------------------
# pollution <- read_csv2("~/Desktop/marine_sediment/data/pollution.csv") %>%
#   dplyr::select(.,-c(`[Water Control Zone]`,STATION,`[Sampling Date]`,`[Sample Number]`,`[Sample Type]`)) %>%
#   group_by(StationID,Month_Ret) %>%
#   summarise_if(is.double,mean) %>%
#   set_names(c("Site","Season","TVS","COD","TC","NH4-N","TKN","Cu","Pb","Zn"))
pollution_file <- read_tsv("~/Desktop/marine_sediment/results/pollution/year2019_pollution.tsv")
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
## OTU table: taxon ------------------------------------------------------------
phylum_file <- paste0(data_path,"phylum.otu")
class_file <- paste0(data_path,"class.otu")
order_file <- paste0(data_path,"order.otu")
family_file <- paste0(data_path,"family.otu")
genus_file <- paste0(data_path,"genus.otu")
species_file <- paste0(data_path,"species.otu")
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
  return(taxon)
}
phylum_otu <- read_taxon(phylum_file)
class_otu <- read_taxon(class_file)
order_otu <- read_taxon(order_file)
family_otu <- read_taxon(family_file)
genus_otu <- read_taxon(genus_file)
species_otu <- read_taxon(species_file)
## Partial Spearman correlation ------------------------------------------------
Partial_Spearman_adjust_season <- function(data,sign)
  {
# data <- phylum_otu
# sign <- ancom2_gradient_list$Phylum$re_1$df$taxa_id

pcor_test <- function(df,x,y,p){
  a <- df[[x]]
  b <- df[[y]]
  c <- df[[p]] %>% as.factor() %>% as.integer()
  r <- as_tibble(pcor.test(a,b,c,method = "spearman"))
  return(r)
}

microbe <- data %>%
  dplyr::select(SampleID,sign) %>%
  separate(SampleID, into = c("Site", "Season","ID"), sep = "_") %>%
  dplyr::select(.,-ID) %>%
  group_by(Site,Season) %>%
  summarise_if(is.double,mean) %>%
  ungroup() %>%
  pivot_longer(3:length(.),names_to = "microbiome",values_to = "abundance") %>%
  mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
  group_by(Site,Season,microbiome) %>%
  summarise(abundance = mean(abundance)) %>%
  ungroup() %>%
  left_join(pollution) %>%
  pivot_longer(5:length(.),names_to = "env_param",values_to = "env_value") %>%
  #filter(metabolite == "Isocitric acid",parameter == "MAP") %>%
  #group_split(metabolite,parameter)
  group_by(microbiome,env_param) %>%
  #do(tidy(pcor.test(After_subtraction,value,`Age(years)`, data = .))
  group_modify(~ pcor_test(.x,"abundance","env_value","Season")) %>%
  bind_rows() %>%
  ungroup() %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) #%>% # p value correction

#corr <- rcorr(as.matrix(Met_mat), as.matrix(Clinical_mat), type = 'spearman')
r <- microbe %>%
  dplyr::select(microbiome,env_param,estimate) %>%
  pivot_wider(names_from = env_param,values_from = estimate) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p <- microbe %>%
  dplyr::select(microbiome,env_param,p.value) %>%
  pivot_wider(names_from = env_param,values_from = p.value) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p.BH <- microbe %>%
  dplyr::select(microbiome,env_param,q.value) %>%
  pivot_wider(names_from = env_param,values_from = q.value) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p_1 <- p
p_1[p >= 0.05] <- 1
p_1[p < 0.05] <- 0
p.BH <- p.BH + p_1
## Heatmap ---------------------------------------------------------------------
re_p <- Heatmap(r,name = "Spearman r",rect_gp = gpar(col = "white", lwd = 2),
        row_names_max_width = max_text_width(rownames(r),gp = gpar(fontsize = 12)),
        column_names_max_height = max_text_width(colnames(r),gp = gpar(fontsize = 12)),
        # label the circle of fdr 
        #layer_fun = function(j, i, x, y, width, height, fill) {
        #  ind_mat = restore_matrix(j, i, x, y)
        #  v = pindex(p.BH, i, j) # qvalue is the matrix : qvalue table 
        #  l = v < 0.2
        #  ind_mat_vec = as.matrix(ind_mat)
        #  ind = ind_mat_vec[l]
        #  grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
        #},
        #right_annotation = bar,
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
          grid.points(x[ind], y[ind], pch = 1, size = unit(2, "mm"), gp = gpar(col = "white",lwd=2))
          
          
          v2 = pindex(p.BH, i, j)
          l2 = v2 < 0.2
          ind_mat_vec2 = as.matrix(ind_mat)
          ind2 = ind_mat_vec2[l2]
          grid.points(x[ind2], y[ind2], pch = 2, size = unit(4, "mm"), gp = gpar(col = "white",lwd=2))
          
        })
  return_list <- list(df = microbe,
                      p = re_p)
  return(return_list)
}
## Sign table ------------------------------------------------------------------
#load("~/Desktop/marine_sediment/RDS_files/ancom2_gradient_list_alpha_0.1_detected_0.7.RData") # ancom2_gradient_list
ancom2_gradient_list <- readRDS("~/Desktop/marine_sediment/RDS_files/ancom2_gradient_list.rds")
Phylum_list <- Partial_Spearman_adjust_season(phylum_otu,ancom2_gradient_list$Phylum$re_1$df %>%
                                                filter(!W %in% "Inf") %>%
                                                filter(detected_0.7 %in% "TRUE") %>% 
                                                .$taxa_id)
#Class_list <- Partial_Spearman_adjust_season(class_otu,ancom2_gradient_list$Class$re_1$df$taxa_id)
Genus_list <- Partial_Spearman_adjust_season(genus_otu,ancom2_gradient_list$Genus$re_1$df %>% 
                                               filter(!W %in% "Inf") %>%
                                               filter(detected_0.7 %in% "TRUE") %>% 
                                               .$taxa_id)
Species_list <- Partial_Spearman_adjust_season(species_otu,ancom2_gradient_list$Species$re_1$df %>%
                                                 filter(!W %in% "Inf") %>%
                                                 filter(detected_0.7 %in% "TRUE") %>% 
                                                 .$taxa_id)
sign_genus <- ancom2_gradient_list$Genus$re_1$df %>% 
  filter(!W %in% "Inf") %>%
  filter(detected_0.7 %in% "TRUE") %>% 
  .$taxa_id
sign_species <- ancom2_gradient_list$Species$re_1$df %>% 
  filter(!W %in% "Inf") %>%
  filter(detected_0.7 %in% "TRUE") %>% 
  .$taxa_id
taxon_info <- read_tsv("~/Desktop/marine_sediment/data/sp2fullinfo.txt")
taxon_tmp <- taxon_info %>%
  separate(Info,c("k", "p","c","o","f","g","s"),sep= ";") %>%
  mutate(g = str_remove_all(g,"g_")) %>%
  filter(g %in% sign_genus) %>%
  distinct(., g, .keep_all = TRUE)
taxon_sign <- taxon_info %>%
  filter(ID %in% sign_species) %>%
  separate(Info,c("k", "p","c","o","f","g","s"),sep= ";")
View(taxon_sign)
# ancom2_season_mle4_list <- readRDS("~/Desktop/marine_sediment/RDS_files/ancom2_season_mle4_list.rds")
# sign_genus <- ancom2_season_mle4_list$Genus$re_1$df %>% 
#   filter(detected_0.7 %in% "TRUE") %>% 
#   .$taxa_id
# write_tsv(Species_list$df,"~/Desktop/marine_sediment/results/Species_env_partial_spearman_correlation.tsv")
sign_tab <- rbind(data.frame(Phylum_list$df,level="Phylum"),
                  data.frame(Genus_list$df,level="Genus"),
                  data.frame(Species_list$df,level="Species")
                  )

r <- sign_tab %>%
  dplyr::select(microbiome,env_param,estimate) %>%
  pivot_wider(names_from = env_param,values_from = estimate) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p <- sign_tab %>%
  dplyr::select(microbiome,env_param,p.value) %>%
  pivot_wider(names_from = env_param,values_from = p.value) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p.BH <- sign_tab %>%
  dplyr::select(microbiome,env_param,q.value) %>%
  pivot_wider(names_from = env_param,values_from = q.value) %>%
  column_to_rownames(var = "microbiome") %>%
  as.matrix()
p_1 <- p
p_1[p >= 0.05] <- 1
p_1[p < 0.05] <- 0
p.BH <- p.BH + p_1

## Heatmap ---------------------------------------------------------------------
Pollution_group <- read_tsv("~/Desktop/marine_sediment/results/pollution/pollution_class.tsv")
Pollution_group$name <- str_remove_all(Pollution_group$name," \\(μg/kg\\)") %>%
  str_remove_all(.," \\(mg/kg\\)") %>%
  str_remove_all(.," \\(mg/L\\)") %>%
  str_remove_all(.," \\(μg/L\\)")
Pollution_group <- Pollution_group %>%
  filter(name %in% colnames(r)) %>%
  mutate(name = factor(name,levels = colnames(r))) %>%
  arrange(name) %>%
  column_to_rownames("name")
all_sign_plot <- sign_tab %>% ungroup() %>%  filter(env_param == "Lead") %>%
  mutate(level = factor(level,levels = c("Phylum","Class","Order","Family","Genus","Species")))
rownames(all_sign_plot) <- all_sign_plot$microbiome
Sign = data.frame(level = all_sign_plot$level)
#cbPalette <- c("Phylum" = "#E69F00", "Class" = "#56B4E9","Order" ="#009E73","Family" = "#CC79A7","Genus"= "#0072B2","Species" = "#D55E00")
cbPalette <- c("Phylum" = "#E69F00","Genus"= "#0072B2","Species" = "#D55E00")
bar <- columnAnnotation(df = Sign,col = list(level =cbPalette),
                        annotation_name_gp = gpar(fontsize = 16),
                     annotation_legend_param = list(nrow = 1,
                                                    labels_gp = gpar(fontsize = 16),
                                                    title_gp = gpar(fontsize = 16,fontface = "bold"))
                     )
pollutionPalette <- c("Heavy metal" ="#555555","Nitrogen" = "#F3C57B","PAH" = "#9A8A76","Organic pollution"="#DB735C")
foo <- rowAnnotation(df = Pollution_group,col = list(group = pollutionPalette),
                     annotation_name_gp = gpar(fontsize = 16),
                        annotation_legend_param = list(title = "Pollution group",
                                                       nrow = 1,
                                                       labels_gp = gpar(fontsize = 16),
                                                       title_gp = gpar(fontsize = 16,fontface = "bold")))
# pdf("~/Desktop/marine_sediment/results/0_background_modify/Sign_correlation.pdf",
#     height = 12,width = 7)
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
pdf("~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiome_v3/Figures_v3/Sign_correlation1.pdf",
    height = 10,width = 23)
draw(ht_list,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()
return(ht_list)
}
Partial_correlation_with_env_from_ANCOM2_list <- Partial_correlation_with_env_from_ANCOM2()
saveRDS(Partial_correlation_with_env_from_ANCOM2_list,
     file="~/Desktop/marine_sediment/RDS_files/Partial_correlation_with_env_from_ANCOM2_list.rds")

