## Load library ----------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
data_path <- "~/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
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
Score_sediment <- df_tidy_sediment %>%
  filter(year %in% c("2019"),season %in% c("Summer","Winter")) %>%
  type_convert() %>%
  group_by(StationID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  dplyr::select(.,-c("year","month","Sample No")) %>%
  #column_to_rownames("StationID") %>%
  dplyr::rename(`TotalPCBs (μg/kg)` = `Total Polychlorinated Biphenyls (μg/kg)`) %>%
  dplyr::select(.,-c(`TotalPCBs (μg/kg)`,`Silver (mg/kg)`))
## load species df -------------------------------------------------------------
#species_df <- read_tsv(paste0(data_path,"species.otu"))
sign_dbRDA_lst <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/dbRDA.rds")
## Heatmap ---------------------------------------------------------------------
sediment_all_names <- read_tsv("~/Desktop/marine_sediment/renew_pollutionbin/sediment_pollute.tsv")
Metal_addINFO <- data.frame(
  name = c("Mercury (mg/kg)","Cadmium (mg/kg)","Arsenic (mg/kg)",
           "Copper (mg/kg)","Lead (mg/kg)","Chromium (mg/kg)",
           "Zinc (mg/kg)","Nickel (mg/kg)"),
  toxicity_response_coefficient = c(40,30,10,5,5,2,2,5),
  weight = c(NA_real_,0.25,NA_real_,0.075,0.251,0.134, 0.075,0.215),
  ERL = c(0.15,5,33,70,35,80,120,30),
  ERM = c(1.3,9,85,390,110,145,270,50)
)

Score_sediment_new <- Score_sediment %>%
  pivot_longer(-1) %>%
  left_join(sediment_all_names)

sign_uni <- sign_dbRDA_lst$R_square_table_cate %>%
  filter(unip <= 0.05, BHP <= 0.05) %>%
  dplyr::filter(!name %in% c("Total Solid (%w/w)","Dry Wet Ratio"))

All_uni_sign <- Score_sediment %>%
  dplyr::select(StationID,sign_uni$name) %>%
  pivot_longer(-1) %>%
  mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
  mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
  pivot_wider() %>%
  column_to_rownames("StationID") %>%
  as.matrix() %>%
  scale() %>%
  t()
Univariate_p <- ComplexHeatmap::Heatmap(All_uni_sign,name = "Normalized value",
        column_names_max_height = max_text_width(
          colnames(All_uni_sign), 
          gp = gpar(fontsize = 16, fontface = "bold")
        ),
        row_names_max_width = max_text_width(
          rownames(All_uni_sign), 
          gp = gpar(fontsize = 16, fontface = "bold")
        ),
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        row_names_gp = gpar(fontsize = 16, fontface = "bold"),
        heatmap_legend_param = list(
          at = c(-3, 0, 3),
          labels = c("low", "zero", "high"),
          title = "Z-scaled values",
          legend_height = unit(4, "cm"),
          title_gp = gpar(fontsize = 16, fontface = "bold"),  # Adjust the legend title font size
          labels_gp = gpar(fontsize = 16, fontface = "bold"),
          title_position = "lefttop-rot"
        )
        )

pdf("~/Desktop/marine_sediment/results/MEGAN6/dbRDA_Univariate_selecrPollutant_siteonly_pollutant.pdf",
    width = 9,height = 5)
Univariate_p
dev.off()

Stepwise_sign <- Score_sediment %>%
  dplyr::select(StationID,sign_dbRDA_lst$forward_select)  %>%
  pivot_longer(-1) %>%
  mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
  mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
  mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
  mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
  pivot_wider() %>%
  column_to_rownames("StationID") %>%
  as.matrix() %>%
  scale() %>%
  t()
Stepwise_p <- ComplexHeatmap::Heatmap(Stepwise_sign,name = "stepwise",
                        column_names_max_height = max_text_width(
                          colnames(Stepwise_sign), 
                          gp = gpar(fontsize = 12)
                        ),
                        row_names_max_width = max_text_width(
                          rownames(Stepwise_sign), 
                          gp = gpar(fontsize = 12)
                        ),
                        heatmap_legend_param = list(
                          at = c(-3, 0, 3),
                          labels = c("low", "zero", "high"),
                          title = "Z-scaled values",
                          legend_height = unit(4, "cm"),
                          title_position = "lefttop-rot"
                        ),
                        clustering_method_columns = "ward.D",
                        clustering_distance_columns  = "euclidean"
)
Stepwise_p
Pollution_df <- as.data.frame(colSums(All_uni_sign)) %>%
  `colnames<-`("Univariate") %>%
  rownames_to_column("Site") %>%
  left_join(as.data.frame(colSums(Stepwise_sign)) %>%
              `colnames<-`("Stepwise") %>%
              rownames_to_column("Site"))

saveRDS(Pollution_df,"~/Desktop/marine_sediment/results/MEGAN6/Pollution_tab.rds") 
  
Pollution_df <- Pollution_df %>%
  mutate(Univariate = round(Univariate, 2)) %>%
  mutate(Stepwise = round(Stepwise, 2))
  
write_tsv(Pollution_df,"~/Desktop/marine_sediment/results/MEGAN6/Pollution_tab_round.tsv")



pdf("~/Desktop/marine_sediment/results/MEGAN6/dbRDA_stepwise_selecrPollutant_siteonly_1.pdf",
    width = 8,height = 5)
Stepwise_p
dev.off()