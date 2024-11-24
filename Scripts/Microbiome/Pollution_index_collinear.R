# Title: Collinear correlation
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

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
## read station sites ----------------------------------------------------------
pollution_sediment <- read_csv2("~/Desktop/marine_sediment/data/pollution.csv") %>%
  select(STATION,StationID) %>%
  distinct() %>%
  dplyr::rename(Station = STATION)
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
remove_pollutants <- remove_pollutants %>%
  filter(name %in% c("Electrochemical Potential (mV)"))
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
  unite("group", StationID:season) %>%
  column_to_rownames("group") %>%
  dplyr::select(.,-`Acenaphthene (μg/kg)`,-`Acenaphthylene (μg/kg)`,
                -`Anthracene (μg/kg)`,-`Fluorene (μg/kg)`) %>%
  as.matrix() 


combo_info <- caret::findLinearCombos(Score_sediment)

library(Hmisc)
corr <- Hmisc::rcorr(Score_sediment)

BH_P <- corr$P %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  pivot_longer(-1) %>%
  mutate(value = replace_na(value,1)) %>%
  mutate(q_value = p.adjust(value,method = "BH")) %>%
  dplyr::select(.,-value) %>%
  pivot_wider(names_from = name, values_from = q_value) %>%
  column_to_rownames("ID")
all(colnames(BH_P) == colnames(corr$P))
all(rownames(BH_P) == rownames(corr$P))
# Define color mapping
col_fun <- circlize::colorRamp2(c(-1,0,1), c("blue","white","red"))

colnames(corr$r) <- str_remove_all(colnames(corr$r)," \\(.*")
rownames(corr$r) <- str_remove_all(rownames(corr$r)," \\(.*")
colnames(corr$P) <- str_remove_all(colnames(corr$P)," \\(.*")
rownames(corr$P) <- str_remove_all(rownames(corr$P)," \\(.*")
colnames(BH_P) <- str_remove_all(colnames(BH_P)," \\(.*")
rownames(BH_P) <- str_remove_all(rownames(BH_P)," \\(.*")

Cor_Ht <- Heatmap(corr$r, name = "Pollution Correlation",
                  col = col_fun, #rect_gp = gpar(type = "none"),
                  rect_gp = gpar(color = "black",fill = NA), 
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    # grid.rect(x = x, y = y, w = w, h = h,
                    #           gp = gpar(col = "black", fill = NA))
                    # grid.circle(x = x, y = y, r = abs(stool_p_mat[i, j])/50 * min(unit.c(w, h)),
                    #             gp = gpar(col = "black", fill=NA))
                    if(corr$r[i,j] >= 0.9) {
                      grid.rect(x, y, w, h, 
                                gp = gpar(fill = col_fun(corr$r[i, j])))
                    }
                    if(BH_P[i,j] <= 0.05 & corr$P[i,j] <= 0.05 & corr$r[i,j] >= 0.9) {
                      grid.points(x , y , pch = 16, size = unit(2, "mm"), gp = gpar(col = "black",lwd=2))
                    }
                  }, 
                  row_names_max_width = max_text_width(
                    rownames(corr$r), 
                    gp = gpar(fontsize = 12)
                  ),
                  column_names_max_height = max_text_width(
                    colnames(corr$r), 
                    gp = gpar(fontsize = 12)
                  ),
                  cluster_rows = T, cluster_columns = T,
                  show_row_names = T, show_column_names = T,
                  column_names_side = "bottom",
                  show_heatmap_legend = T
)
Cor_Ht

cairo_pdf("~/Desktop/marine_sediment/results/Pollution_indices_correlation.pdf",
    width = 12,height = 12)
Cor_Ht
dev.off()




