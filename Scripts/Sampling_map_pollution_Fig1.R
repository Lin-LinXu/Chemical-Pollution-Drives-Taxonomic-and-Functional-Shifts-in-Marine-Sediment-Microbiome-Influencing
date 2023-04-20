## Load library ----------------------------------------------------------------
library(tidyverse)
library(Hmisc) ## calculate microbiom correlation coefficient
library(psych)
library(spearmanCI)
library(vegan)
library('RColorBrewer')
library(ComplexHeatmap)
data_path <- "~/Desktop/marine_sediment/data/"
## load sediment data -------------------------------------------------------------------
file1 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ms5.csv")
file2 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ns6.csv")
file3 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ps6.csv")
file4 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ss1-ss3-ss5.csv")
file5 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ts2.csv")
file6 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5yrs-ms8.csv")
file_sediment <- rbind(file1,file2,file3,file4,file5,file6)
## load water data -------------------------------------------------------------------
file1 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_5y_mm5-mm8.csv")
file2 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_5y_nm8.csv")
file3 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_5y_pm11_pm9.csv")
file4 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_5y_sm1-sm18-sm5-sm19-sm9-sm10.csv")
file5 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_5y_tm4.csv")
file_water <- rbind(file1,file2,file3,file4,file5)
Water_depth <- read_tsv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_water_quality_station_depth.tsv")
file_water_tmp <- Water_depth %>%
  left_join(file_water)
## pollution data: pollution ---------------------------------------------------
pollution_sediment <- read_csv2("~/Desktop/marine_sediment/data/pollution.csv") %>%
  select(STATION,StationID) %>%
  distinct() %>%
  rename(Station = STATION)
pollution_water <- read_tsv("~/Desktop/marine_sediment/2019 Water and Sediment Data/MARINE_WATER_DATA.tsv") %>%
  select(Station,StationID) %>%
  distinct() %>%
  mutate(StationID = str_replace(StationID,"AIR","SSW"))
df_sediment <- file_sediment %>%
  left_join(pollution_sediment) %>%
  relocate(StationID,.before = Station)
df_water <- file_water_tmp
gradient <- tibble(StationID = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                   Gradient = c("High","High","High","Medium","Medium","Low","Low")) %>%
  mutate(StationID = factor(StationID,levels = c("SSW","PC","CI","SW","CDA","BI","TPC")))
df_tidy_sediment <- df_sediment %>%
  # select(.,-c("Total Cyanide (mg/kg)","Silver (mg/kg)","Naphthalene (μg/kg)",
  #             "Fluorene (μg/kg)","Dibenzo(ah)anthracene (μg/kg)",
  #             "Chrysene (μg/kg)","Acenaphthene (μg/kg)","Acenaphthylene (μg/kg)",
  #             "Anthracene (μg/kg)")) %>%
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
df_tidy_water <- df_water %>%
  mutate_all(funs(str_replace(., "<", ""))) %>%
  mutate_all(funs(str_replace(., "N/A", NA_character_))) %>%
  as_tibble() %>%
  type_convert() %>%
  separate(Dates,c("year","month","day"), sep = "-") %>%
  mutate(month = as.numeric(month)) %>%
  mutate(
    season = case_when(
      month %in% 10:12 ~ "Fall",
      month %in%  1:3  ~ "Winter",
      month %in%  4:6  ~ "Spring",
      TRUE ~ "Summer")) %>%
  #filter(season %in% c("Summer","Winter")) %>%
  relocate(season,.before = month) #%>%
  #mutate(DIN = `Ammonia Nitrogen (mg/L)` + `Nitrite Nitrogen (mg/L)` + `Nitrate Nitrogen (mg/L)`)
## calculate score -------------------------------------------------------------
Score_sediment <- df_tidy_sediment %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Type")) %>%
  dplyr::select("StationID","Copper (mg/kg)","Dibenzo(ah)anthracene (μg/kg)","Indeno(1,2,3-cd)pyrene (μg/kg)",
                "Benzo(ghi)perylene (μg/kg)","Fluoranthene (μg/kg)","Benzo(a)anthracene (μg/kg)","Benzo(a)pyrene (μg/kg)",
                "Mercury (mg/kg)","Benzo(b)fluoranthene (μg/kg)","Benzo(k)fluoranthene (μg/kg)",
                "Ammonia Nitrogen (mg/kg)","Zinc (mg/kg)","Lead (mg/kg)","Cadmium (mg/kg)",
                #"Total Sulphide (mg/kg)",
                "Total Kjeldahl Nitrogen (mg/kg)","Arsenic (mg/kg)",
                "Total Volatile Solid (% Total Solid)","Chemical Oxygen Demand (mg/kg)",
                "Total Carbon (%w/w)","Total Phosphorus (mg/kg)",
                "Manganese (mg/kg)","Aluminium (mg/kg)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  pivot_longer(2:length(.)) %>%
  group_by(name) %>%
  mutate(value = if_else(value >= quantile(value,probs = 2/3,na.rm =T),3
                         ,if_else(value < quantile(value,probs = 1/3,na.rm =T),1,2))) %>%
  ungroup() %>%
  group_by(StationID) %>%
  summarise(value = sum(value)) %>%
  mutate(level = if_else(value >= quantile(value,probs = 2/3),3
                         ,if_else(value < quantile(value,probs = 1/3),1,2)))

Score_water <- df_tidy_water %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Depth")) %>%
  rename(StationID = "Site") %>%
  dplyr::select("StationID","Turbidity (NTU)","Suspended Solids (mg/L)","Silica (mg/L)","Nitrate Nitrogen (mg/L)",
                "Total Inorganic Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Volatile Suspended Solids (mg/L)",
                "Phaeo-pigments (μg/L)",
                #"Chlorophyll-a (μg/L)",
                "Orthophosphate Phosphorus (mg/L)",
                "Total Nitrogen (mg/L)","Total Phosphorus (mg/L)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  pivot_longer(2:length(.)) %>%
  group_by(name) %>%
  mutate(value = if_else(value >= quantile(value,probs = 2/3),3
                         ,if_else(value < quantile(value,probs = 1/3),1,2))) %>%
  ungroup() %>%
  group_by(StationID) %>%
  summarise(value = sum(value)) %>%
  mutate(level = if_else(value >= quantile(value,probs = 2/3),3
                         ,if_else(value < quantile(value,probs = 1/3),1,2)))

Final_score <- data.frame(
  name = Score_sediment$StationID,
  score = Score_sediment$value + Score_water$value
) %>%
  mutate(level = if_else(score >= quantile(score,probs = 2/3),3
                         ,if_else(score < quantile(score,probs = 1/3),1,2)))

Score_sediment <- df_tidy_sediment %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Type")) %>%
  dplyr::select("StationID","Copper (mg/kg)","Dibenzo(ah)anthracene (μg/kg)","Indeno(1,2,3-cd)pyrene (μg/kg)",
                "Benzo(ghi)perylene (μg/kg)","Fluoranthene (μg/kg)","Benzo(a)anthracene (μg/kg)","Benzo(a)pyrene (μg/kg)",
                "Mercury (mg/kg)","Benzo(b)fluoranthene (μg/kg)","Benzo(k)fluoranthene (μg/kg)",
                "Ammonia Nitrogen (mg/kg)","Zinc (mg/kg)","Lead (mg/kg)","Cadmium (mg/kg)",
                #"Total Sulphide (mg/kg)",
                "Total Kjeldahl Nitrogen (mg/kg)","Arsenic (mg/kg)",
                "Total Volatile Solid (% Total Solid)","Chemical Oxygen Demand (mg/kg)",
                "Total Carbon (%w/w)","Total Phosphorus (mg/kg)",
                "Manganese (mg/kg)","Aluminium (mg/kg)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  pivot_longer(2:length(.)) %>%
  group_by(name) %>%
  mutate(value = if_else(value >= quantile(value,probs = 2/3,na.rm =T),3
                         ,if_else(value < quantile(value,probs = 1/3,na.rm =T),1,2))) %>%
  ungroup()
Score_water <- df_tidy_water %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Depth")) %>%
  rename(StationID = "Site") %>%
  dplyr::select("StationID","Turbidity (NTU)","Suspended Solids (mg/L)","Silica (mg/L)","Nitrate Nitrogen (mg/L)",
                "Total Inorganic Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Volatile Suspended Solids (mg/L)",
                "Phaeo-pigments (μg/L)",
                #"Chlorophyll-a (μg/L)",
                "Orthophosphate Phosphorus (mg/L)",
                "Total Nitrogen (mg/L)","Total Phosphorus (mg/L)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  pivot_longer(2:length(.)) %>%
  group_by(name) %>%
  mutate(value = if_else(value >= quantile(value,probs = 2/3),3
                         ,if_else(value < quantile(value,probs = 1/3),1,2))) %>%
  ungroup()
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
Score_sites <- rbind(Score_sediment,Score_water) %>%
  ungroup() %>%
  count(StationID, value, sort = TRUE) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(percentage = (n/31) )
Label_df <- Score_sites %>%
  dplyr::select(StationID, value,percentage) %>%
  pivot_wider(names_from = "value",values_from = "percentage") %>%
  `colnames<-`(c("StationID","High","Low","Medium")) %>%
  mutate(Low = replace_na(Low,0))
map_id <- data.frame(
  StationID = c("SSW","PC","CI","SW","CDA","BI","TPC"),
  long = c(113.898467,114.017967,114.219600,114.079050,114.231417,114.312417,114.393883),
  lat = c(22.339433,22.257383,22.432733,22.191667,22.212300,22.319467,22.520550)
) %>%
  mutate(StationID = factor(StationID,levels = StationID))
Label_df <- Label_df %>%
  left_join(map_id)
Label_df$long <- Label_df$long + 0.02
library(ggmap)
library(scatterpie)
register_google(key = "AIzaSyAfqtVWkUTf8f4itKHZ0-WBObIp7gKZ0AY")
map <- get_googlemap(center = 'Hong Kong',zoom = 10, scale = 2, maptype = "terrain", 
                     style = 'feature:all|element:labels|visibility:off')
saveRDS(map,"~/Desktop/marine_sediment/RDS_files/map.rds")
cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
Final_map <- ggmap(map) +
  labs(x= "",y="") +
  # geom_point(aes(x = long, y = lat,color = StationID), 
  #            data = map_id, size=2) +
  # geom_text(aes(x = long, y = lat,label = StationID,color = StationID),data = map_id,  
  #           position = position_dodge(width = 0),
  #           vjust = -0.5,size = 4) +
  scale_colour_manual(values = cbPalette) + 
  geom_scatterpie(aes(x=long, y=lat, group=StationID),
                  pie_scale = 2.5,
                    data=Label_df, cols=c("Low","Medium","High"),color = NA, alpha=1) +
  scale_fill_manual(values=c("#0072B2","#E69F00","#D55E00")) +
  theme(legend.position="none")
pdf("~/Desktop/marine_sediment/manuscript/Figure_v5//Map.pdf",
    width = 12,height = 10)
Final_map
dev.off()
#write_csv(Score_sediment,"~/Desktop/marine_sediment/results/pollution/Score_sediment.csv")
#write_csv(Score_water,"~/Desktop/marine_sediment/results/pollution/Score_water.csv")
#write_csv(Final_score,"~/Desktop/marine_sediment/results/pollution/Score_sediment_and_water.csv")
##  -------------------------------------------------------------

# tmp <- df_tidy_sediment %>%
#   filter(year == "2019") %>%
#   dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Type")) %>%
#   group_by(StationID) %>%
#   summarise_if(is.numeric, mean, na.rm = TRUE) %>%
#   column_to_rownames("StationID") %>%
#   caret::nearZeroVar()
Tab_S1 <- df_tidy_sediment %>%
  dplyr::select("StationID","Water Control Zone","Station","year","season","month","day","Sample No","Type",
                "Copper (mg/kg)","Dibenzo(ah)anthracene (μg/kg)","Indeno(1,2,3-cd)pyrene (μg/kg)",
                "Benzo(ghi)perylene (μg/kg)","Fluoranthene (μg/kg)","Benzo(a)anthracene (μg/kg)","Benzo(a)pyrene (μg/kg)",
                "Mercury (mg/kg)","Benzo(b)fluoranthene (μg/kg)","Benzo(k)fluoranthene (μg/kg)",
                "Ammonia Nitrogen (mg/kg)","Zinc (mg/kg)","Lead (mg/kg)","Cadmium (mg/kg)",
                #"Total Sulphide (mg/kg)",
                "Total Kjeldahl Nitrogen (mg/kg)","Arsenic (mg/kg)",
                "Total Volatile Solid (% Total Solid)","Chemical Oxygen Demand (mg/kg)",
                "Total Carbon (%w/w)","Total Phosphorus (mg/kg)",
                "Manganese (mg/kg)","Aluminium (mg/kg)")
Tab_S2 <- df_tidy_water %>%
  rename(StationID = "Site") %>%
  dplyr::select("StationID","Water Control Zone","Station","year","season","month","day","Sample No","Depth",
                "Turbidity (NTU)","Suspended Solids (mg/L)","Silica (mg/L)","Nitrate Nitrogen (mg/L)",
                "Total Inorganic Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Volatile Suspended Solids (mg/L)",
                "Phaeo-pigments (μg/L)",
                #"Chlorophyll-a (μg/L)",
                "Orthophosphate Phosphorus (mg/L)",
                "Total Nitrogen (mg/L)","Total Phosphorus (mg/L)")
write_tsv(Tab_S1,"~/Desktop/marine_sediment/manuscript/Figure_v5/Tab_S1.tsv")
write_tsv(Tab_S2,"~/Desktop/marine_sediment/manuscript/Figure_v5/Tab_S2.tsv")

Mat_sediment <- df_tidy_sediment %>%
  #filter(year == "2019") %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Type")) %>%
  dplyr::select("StationID","Copper (mg/kg)","Dibenzo(ah)anthracene (μg/kg)","Indeno(1,2,3-cd)pyrene (μg/kg)",
                "Benzo(ghi)perylene (μg/kg)","Fluoranthene (μg/kg)","Benzo(a)anthracene (μg/kg)","Benzo(a)pyrene (μg/kg)",
                "Mercury (mg/kg)","Benzo(b)fluoranthene (μg/kg)","Benzo(k)fluoranthene (μg/kg)",
                "Ammonia Nitrogen (mg/kg)","Zinc (mg/kg)","Lead (mg/kg)","Cadmium (mg/kg)",
                #"Total Sulphide (mg/kg)",
                "Total Kjeldahl Nitrogen (mg/kg)","Arsenic (mg/kg)",
                "Total Volatile Solid (% Total Solid)","Chemical Oxygen Demand (mg/kg)",
                "Total Carbon (%w/w)","Total Phosphorus (mg/kg)",
                "Manganese (mg/kg)","Aluminium (mg/kg)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  #dplyr::select(.,-"Acenaphthylene (μg/kg)") %>%
  # pivot_longer(2:length(.)) %>%
  # group_by(name) %>%
  # mutate(value = scale(value)) %>%
  # pivot_wider() %>%
  column_to_rownames("StationID") %>%
  scale() %>%
  as.matrix()
#Mat_sediment <- Mat_sediment[,-tmp]
Mat_water <- df_tidy_water %>%
  dplyr::select(.,-c("Water Control Zone","Station","year","season","month","day","Sample No","Depth")) %>%
  rename(StationID = "Site") %>%
  dplyr::select("StationID","Turbidity (NTU)","Suspended Solids (mg/L)","Silica (mg/L)","Nitrate Nitrogen (mg/L)",
                "Total Inorganic Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Volatile Suspended Solids (mg/L)",
                "Phaeo-pigments (μg/L)",
                #"Chlorophyll-a (μg/L)",
                "Orthophosphate Phosphorus (mg/L)",
                "Total Nitrogen (mg/L)","Total Phosphorus (mg/L)") %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  column_to_rownames("StationID") %>%
  scale()
Mat_all <- cbind(Mat_sediment,Mat_water) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(rowname = factor(rowname,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
  arrange(rowname) %>%
  column_to_rownames() %>%
  as.matrix() #%>%
  #t()
Mat_plot <- t(Mat_all)
rownames(Mat_plot) <- str_remove_all(rownames(Mat_plot),"\\(μg/kg\\)") %>%
  str_remove_all(.,"\\(mg/kg\\)") %>%
  str_remove_all(.,"\\(mg/L\\)") %>%
  str_remove_all(.,"\\(μg/L\\)")
pdf("~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiome_v3/Figures_v3/Mat_all_selected_HC.pdf",
    width = 6,height = 12)
name_df <- data.frame(Environment = c(rep("sediment",dim(Mat_sediment)[2]),
                                rep("water",dim(Mat_water)[2]))) %>%
  `rownames<-`(c(colnames(Mat_sediment),colnames(Mat_water))) %>%
  as.vector()
envPalette <- c("sediment" ="#B82E2E","water" = "#22AA99")
bar <- rowAnnotation(df = name_df,col = list(Environment = envPalette),
                     show_annotation_name = F,
                     annotation_legend_param = list(title = "Environment",
                                                    nrow = 1,labels_gp = gpar(fontsize = 16),
                                                    title_gp = gpar(fontsize = 16,fontface = "bold")))
col_df <- gradient %>%
  column_to_rownames("StationID") %>%
  as.vector() %>%
  `colnames<-`("Pollution gradient")
col_df$`Pollution gradient` <- factor(col_df$`Pollution gradient`,levels = c("High","Medium","Low"))
pollutionPalette <- c("High" ="#D55E00","Medium" = "#E69F00","Low" = "#0072B2")
foo <- columnAnnotation(df = col_df,col = list(`Pollution gradient` = pollutionPalette),
                        annotation_name_gp = gpar(fontsize = 16),
                        annotation_legend_param = list(title = "Pollution gradient",
                                                       nrow = 1,labels_gp = gpar(fontsize = 16),
                                                       title_gp = gpar(fontsize = 16,fontface = "bold")))
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#01befe", "white", "#ff7d00"))
Ht = ComplexHeatmap::Heatmap(Mat_plot,cluster_columns = T,
                        rect_gp = gpar(col = "black", lwd = 2),
                        right_annotation = bar,
                        top_annotation = foo,
                        col = col_fun,
                        column_names_gp = grid::gpar(fontsize = 16),
                        row_names_gp = grid::gpar(fontsize = 16),
                        row_names_max_width = max_text_width(rownames(Mat_plot),gp = gpar(fontsize = 16)),
                        column_names_max_height = max_text_width(colnames(Mat_plot),gp = gpar(fontsize = 16)),
                        heatmap_legend_param = list(title = "Normalized abundance",direction = "horizontal",
                                                    labels_gp = gpar(fontsize = 16),
                                                    title_gp = gpar(fontsize = 16,fontface = "bold")))
draw(Ht,merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()
# pdf("~/Desktop/marine_sediment/results/pollution/Mat_all_selected.pdf",width = 10,height = 12)
# name_df <- data.frame(group = c(rep("sediment",dim(Mat_sediment)[2]),
#                                 rep("water",dim(Mat_water)[2]))) %>%
#   `rownames<-`(c(colnames(Mat_sediment),colnames(Mat_water))) %>%
#   as.vector()
# envPalette <- c("sediment" ="#FF9900","water" = "#3366CC")
# bar <- rowAnnotation(df = name_df,col = list(group = envPalette))
# col_df <- gradient %>%
#   column_to_rownames("StationID") %>%
#   as.vector()
# pollutionPalette <- c("High" ="#D55E00","Medium" = "#E69F00","Low" = "#0072B2")
# foo <- columnAnnotation(df = col_df,col = list(Gradient = pollutionPalette))
# 
# ComplexHeatmap::Heatmap(Mat_plot,cluster_columns = F,
#                         rect_gp = gpar(col = "black", lwd = 2),
#                         right_annotation = bar,
#                         top_annotation = foo,
#                         row_names_max_width = max_text_width(rownames(Mat_plot),gp = gpar(fontsize = 12)),
#                         column_names_max_height = max_text_width(colnames(Mat_plot),gp = gpar(fontsize = 12)))
# dev.off()
## PCA
####Euclidean distance
Euclidean_dist <- vegdist(Mat_all, method = "euclidean")
Euclidean_dist <- as.matrix(Euclidean_dist)
SSw = 0.5 * 2.609442 + 0.5* 4.342632 + 1/3 * (7.504188+ 10.554216+ 8.226012)
SSt = 1/7 * 159.8605
SSa = SSt - SSw
(SSa - SSw/2)/(SSt + SSw/4)
###PERMANOVA
adonis_Euclidean_dist = adonis2(Euclidean_dist~gradient$Gradient,permutations = 9999)
adonis_Euclidean_R2 = round(adonis_Euclidean_dist$R2[1],2)
adonis_Euclidean_P = round(adonis_Euclidean_dist$`Pr(>F)`[1],4)
# ## PERMANOVA High VS medium
# mat_high_medium <- Mat_all[c("SSW","PC","CI","SW","CDA"),]
# Euclidean_dist <- vegdist(mat_high_medium, method = "euclidean")
# Euclidean_dist <- as.matrix(Euclidean_dist)
# ###PERMANOVA
# adonis_Euclidean_dist = adonis2(Euclidean_dist~c("High","High","High","Medium","Medium"),permutations = 9999)
# adonis_Euclidean_R2 = round(adonis_Euclidean_dist$R2[1],2)
# adonis_Euclidean_P = round(adonis_Euclidean_dist$`Pr(>F)`[1],4)
## P = 0.4
# ## PERMANOVA High VS low
# mat_high_low <- Mat_all[c("SSW","PC","CI","BI","TPC"),]
# Euclidean_dist <- vegdist(mat_high_low, method = "euclidean")
# Euclidean_dist <- as.matrix(Euclidean_dist)
# ###PERMANOVA
# adonis_Euclidean_dist = adonis2(Euclidean_dist~c("High","High","High","Low","Low"),permutations = 9999)
# adonis_Euclidean_R2 = round(adonis_Euclidean_dist$R2[1],2)
# adonis_Euclidean_P = round(adonis_Euclidean_dist$`Pr(>F)`[1],4)
## P = 0.2
####################################################
pca.out <- prcomp(Mat_all)
pr.var<-pca.out$sdev^2
pve<-pr.var/sum(pr.var)
title = sprintf("Principle Component Analysis (%s)",name)
pc1 <- paste("PC1 (",round(pve[1]*100,digits=2),"%)",sep="")
pc2 <- paste("PC2 (",round(pve[2]*100,digits=2),"%)",sep="")
pca <- data.frame(PC1=pca.out$x[,1],PC2=pca.out$x[,2],Site=gradient$StationID,
                  Gradient = gradient$Gradient)
#line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Euclidean_R2)~"P = "~.(adonis_Euclidean_P))
#line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Euclidean_P))
#line_3 = bquote("ANOSIM:"~"R"~"= "~.(anosim_Euclidean_R)~"P = "~.(anosim_Euclidean_P))
#nonmetric_label = list(line_1)
#coord_x <- Inf#max(pca$PC1)
#coord_y <- c(max(pca$PC2))
### Setting color parameters ---------------------------------------------------
library(ggprism)
cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
envPalette <- rev(c("#D95F02","#E6AB02","#1B9E77"))
pca_plot <- ggplot(pca, aes(x=PC1,y=PC2))+
  geom_point(aes(color = Site),position = "jitter",size = 4) +
  labs(x=pc1,y=pc2) + 
  #scale_shape_manual(values = c(17,16))+
  scale_colour_manual(values = cbPalette) + 
  #labs(title = title)+
  # annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),
  #          hjust = 1, vjust = 1,size = 6) +
  #theme_bw()+
  theme_prism(base_size = 16)+
  # stat_ellipse(aes(group = Gradient,fill=Gradient),
  #              geom="polygon",level=0.68,alpha=0.1,
  #              type = "norm")+
  scale_fill_manual(values = envPalette)+
  #theme_prism(border = TRUE,base_size = 16,base_rect_size = 2) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position="none")
pca_plot
saveRDS(pca_plot,"~/Desktop/marine_sediment/results_new/pca_pollution.rds")
library(patchwork)
Final_map + inset_element(pca_plot, left = 0.5, bottom = 0.05, right = 0.9, top = 0.45)
pdf("~/Desktop/marine_sediment/manuscript/Figure_v2/Map.pdf",
    width = 12,height = 10)
Final_map
dev.off()
pdf("~/Desktop/marine_sediment/manuscript/Figure_v5/Pollution_PCA.pdf",
    width = 6,height = 4)
pca_plot
dev.off()
