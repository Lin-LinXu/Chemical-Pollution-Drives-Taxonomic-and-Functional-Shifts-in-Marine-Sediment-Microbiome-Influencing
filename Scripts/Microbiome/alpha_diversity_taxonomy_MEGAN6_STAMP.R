# Title: Alpha diversity microbiome
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

### Loading required package ---------------------------------------------------
library(vegan)
library(ggplot2)
library(reshape)
library(dplyr)
library(vegan)
library(ggpubr)
library(patchwork)
library(rstatix)
library(numbers)
library(tidyverse)
library(ggprism)
library(phyloseq)
data_path <- "/data/Desktop/marine_sediment/MEGAN6/STAMP_OTU/"
### Setting color parameters ---------------------------------------------------
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
Rarefy_stool_microbiome_data <- function(dataframe) {
  # dataframe <- taxon
  
  df <- dataframe %>%
    dplyr::rename(Feature = 1) %>%
    column_to_rownames("Feature") %>%
    as.matrix()
  physeq.species <- phyloseq(otu_table(df,taxa_are_rows = TRUE))
  set.seed(123)
  rarefy.species <- rarefy_even_depth(physeq.species)
  OTU_tab <- data.frame(otu_table(rarefy.species)) %>%
    rownames_to_column("Feature") 
  return(OTU_tab)
}
Flt_taxa <- function(dataframe){
  
  dataframe <- read_tsv(dataframe)
  otu_tab <- dataframe %>%
    dplyr::rename(Feature = 1) %>%
    pivot_longer(-1) %>%
    group_by(Feature) %>%
    summarise(NonZeroCount = sum(value != 0)) %>%
    filter(NonZeroCount > 1)
  final_tab <- dataframe %>%
    dplyr::rename(Feature = 1) %>%
    filter(Feature %in% otu_tab$Feature)
  return(final_tab)
}
Alpha_diversity <- function(data,name){
  # data <- paste0(data_path,"species.otu")
  # name <- "Species"
  
  taxon <- read_tsv(data)
  #taxon <- Flt_taxa(data)
  rarefied_taxon <- Rarefy_stool_microbiome_data(taxon) %>%
    filter(!Feature %in% "Unclassified") %>%
    column_to_rownames("Feature") %>%
    t()
  
  richness_df <- estimateR(rarefied_taxon) %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::select(S.obs,S.chao1,S.ACE) %>%
    `colnames<-`(c("Richness","Chao1","ACE")) %>%
    rownames_to_column("SampleID")
  
  taxon1 <- Renormalize_data(taxon) %>%
    filter(!Feature %in% "Unclassified") %>%
    column_to_rownames(var = "Feature") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
  
  shannon_diversity <- as.data.frame(diversity(taxon1)) %>%
    `colnames<-`("Shannon_diversity")
  simpson_diversity <- as.data.frame(diversity(taxon1,index="simpson")) %>%
  `colnames<-`("Simpson_diversity")
  alpha_diversity <- merge(shannon_diversity,simpson_diversity,by="row.names",
                           all.x = TRUE) %>%
    dplyr::rename(SampleID = Row.names) %>%
    left_join(richness_df) %>%
    mutate(tmp = SampleID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble()
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  alpha_diversity <- alpha_diversity %>%
    left_join(gradient)
  
  alpha_diversity_gradient <- alpha_diversity %>%
    dplyr::select(.,-ACE) %>%
    pivot_longer(Shannon_diversity:Chao1) %>%
    group_by(name) %>%
    wilcox_test(value ~ Gradient)
  
  # The palette with grey:
  cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02", "CI" = "#e7298a", 
                 "SW" = "#b10c0c", "CDA" = "#ab37c9", "BI" = "#e96e0f",
                 "TPC" = "#019a5b"
  )
  
  alpha_diversity_season <- alpha_diversity %>%
    dplyr::select(.,-ACE) %>%
    pivot_longer(Shannon_diversity:Chao1) %>%
    group_by(Site,Season,name) %>%
    summarise_at("value",mean) %>%
    ungroup() %>%
    group_by(name) %>%
    wilcox_test(value ~ Season,paired = T)
  
  alpha_diversity_all <- alpha_diversity %>%
    pivot_longer(Shannon_diversity:Chao1) %>%
    mutate(Site = factor(Site, levels = c("SSW", "PC", "CI", 
                                          "SW", "CDA", "BI",
                                          "TPC")
    )
    )
  
  p0 <- ggboxplot(alpha_diversity_all, x = "Site", y = "value",outlier.shape=NA) +
    geom_point(data = alpha_diversity_all, aes(x = Site, y = value, shape=Season,color = Site),
               position = position_jitter(width = 0.2))+
    facet_wrap(.~name,scales = "free") +
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette)
  

  gradient_p <- alpha_diversity_gradient %>% filter(name %in% "Shannon_diversity") %>%
    .$p %>% round(2)
  season_p <- alpha_diversity_season %>% filter(name %in% "Shannon_diversity") %>%
    .$p %>% round(2)
  gradient.label <- paste("Pollution P-value:",gradient_p)
  season.label <- bquote("Season P-value:"~.(season_p))
  nonmetric_label = list(gradient.label,season.label)
  diff_num <- max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Shannon_diversity) + 0.2 * diff_num,
               max(alpha_diversity$Shannon_diversity) + 0.1 * diff_num )
  
  p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Shannon_diversity",outlier.shape=NA) +
    ylab(paste("Shannon diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Shannon_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    ggprism::theme_prism(base_size = 16)
  p1$theme[c("legend.text.align", "legend.title.align")] <- NULL
  #p1
  p2 <- ggboxplot(alpha_diversity, x = "Season", y = "Shannon_diversity",outlier.shape=NA) +
    ylab(paste("Shannon diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Season, y = Shannon_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    theme_prism(base_size = 16)
  p2$theme[c("legend.text.align", "legend.title.align")] <- NULL
  
  gradient_p <- alpha_diversity_gradient %>% filter(name %in% "Simpson_diversity") %>%
    .$p %>% round(2)
  season_p <- alpha_diversity_season %>% filter(name %in% "Simpson_diversity") %>%
    .$p %>% round(2)
  gradient.label <- paste("Pollution P-value:",gradient_p)
  season.label <- bquote("Season P-value:"~.(season_p))
  nonmetric_label = list(gradient.label,season.label)
  diff_num <- max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Simpson_diversity) + 0.2 * diff_num,
               max(alpha_diversity$Simpson_diversity) + 0.1 * diff_num )
  
  p3 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Simpson_diversity",outlier.shape=NA) +
    ylab(paste("Simpson diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Simpson_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    theme_prism(base_size = 16)
  p3$theme[c("legend.text.align", "legend.title.align")] <- NULL
  #p1
  p4 <- ggboxplot(alpha_diversity, x = "Season", y = "Simpson_diversity",outlier.shape=NA) +
    ylab(paste("Simpson diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Season, y = Simpson_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    theme_prism(base_size = 16)
  p4$theme[c("legend.text.align", "legend.title.align")] <- NULL
  
  
  return_result <- list(Shannon=p1,
                        Shannon_season = p2,
                        Simpson=p3,
                        Simpson_season = p4,
                        all_p = p0,
                        sign_gradient = alpha_diversity_gradient,
                        sign_season = alpha_diversity_season,
                        alpha_diversity = alpha_diversity)
  return(return_result)
}
#### Main part -----------------------------------------------------------------
species_list <- Alpha_diversity(paste0(data_path,"species.otu"),
                                "Species"
                                )

genus_list <- Alpha_diversity(paste0(data_path,"genus.otu"),
                                "Genus"
)


saveRDS(species_list,"/data/Desktop/marine_sediment/results/MEGAN6/alpha_diversity_species_STAMP.rds")

saveRDS(genus_list,"/data/Desktop/marine_sediment/results/MEGAN6/alpha_diversity_genus_STAMP.rds")
