# Title: Taxonomic profiling
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

### Loading required package ---------------------------------------------------
library(vegan)
library(ggplot2)
library(reshape)
library(dplyr)
library(patchwork)
library(tidyverse)
library(randomcoloR)
library(janitor)
library(tidyverse)
library(phyloseq)
data_path <- "~/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
### Setting color parameters ---------------------------------------------------
manualcolors <- c('forestgreen','red2','orange','cornflowerblue','magenta',
                  'darkolivegreen4','indianred1','tan4','darkblue','mediumorchid1',
                  'firebrick4','yellowgreen','lightsalmon','tan3','tan1',
                  'darkgray','wheat4','#DDAD4B','chartreuse','seagreen1',
                  'moccasin','mediumvioletred','seagreen','cadeblue1','darkolivegreen1',
                  'tan2','tomato3','#7CE3D8','gainsboro','black')

#### Stack Plot function -------------------------------------------------------
Renormalize_data <- function(dataframe) {
  # dataframe <- Procrustes_sign$Profile1[[1]]
  
  otu_tab <- dataframe %>%
    dplyr::rename(Feature = 1) %>%
    column_to_rownames("Feature")
  GPr  = transform_sample_counts(otu_table(otu_tab,taxa_are_rows = TRUE), function(x) x / sum(x) )
  OTU_tab <- data.frame(otu_table(GPr)) %>%
    rownames_to_column("Feature")
  
  return(OTU_tab)
}   
Stack_plot_season_with_others <- function(data,name){
  # data <- paste0(data_path,"kingdom.otu")
  # name <- "Phylum"
  taxon_tmp <- Renormalize_data(read_tsv(data)) %>%
    dplyr::rename(SampleID = 1) %>%
    dplyr::filter(!SampleID %in% "Unclassified") %>%
    mutate(SampleID = str_remove_all(SampleID,"\\<phylum\\>")) %>%
    mutate(mean = rowMeans(across(where(is.numeric)))) %>%
    arrange(desc(mean)) %>%
    slice_head(n = 20)
  levels_order <- c(taxon_tmp$SampleID,"others")
  taxon <- taxon_tmp %>%
    mutate(SampleID = factor(SampleID)) %>%
    dplyr::select(.,-mean) %>%
    column_to_rownames(var = "SampleID") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleID") %>%
    mutate(tmp = SampleID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble() %>%
    mutate(others = 1 - rowSums(across(where(is.numeric)))) %>%
    relocate(others, .after = SampleID) %>%
    pivot_longer(2:(length(.)-2),names_to = "variable",values_to = "value") %>%
    mutate(variable = factor(variable,levels = levels_order)) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))) %>%
    arrange(Site) 
  taxon <- taxon %>%
    mutate(SampleID = factor(SampleID,levels = unique(taxon$SampleID)))
  ## Using Site, SampleID as id variables
  barchart1 <- ggplot(taxon, aes(x = Site, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Season, scales = 'free_x', ncol = 2) +
    theme_bw() +
    scale_fill_manual(name = name, values = manualcolors) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
          strip.background = element_rect(color="black", 
                                          fill="white",
                                          size= 0.3, 
                                          linetype="solid")
    )
  barchart2 <- ggplot(taxon, aes(x = SampleID, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Season, scales = 'free_x', ncol = 2) +
    theme_bw() +
    scale_fill_manual(name = name, values = manualcolors) +
    labs(x = "", y = "Percentage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          strip.text.x = element_text(size = 12, colour = "black", face = "bold"),
          strip.background = element_rect(color="black", 
                                          fill="white",
                                          size= 0.3, 
                                          linetype="solid")
    )
  return_list <- list(
    barchart_site = barchart1,
    barchart_sample = barchart2
  )
  return(return_list)
}
#### plot taxonomic profiling --------------------------------------------------
phylum_list <- Stack_plot_season_with_others(paste0(data_path,"phylum.otu"),"Phylum")
family_list <- Stack_plot_season_with_others(paste0(data_path,"family.otu"),"Family")
genus_list <- Stack_plot_season_with_others(paste0(data_path,"genus.otu"),"Genus")
species_list <- Stack_plot_season_with_others(paste0(data_path,"species.otu"),"Species")


