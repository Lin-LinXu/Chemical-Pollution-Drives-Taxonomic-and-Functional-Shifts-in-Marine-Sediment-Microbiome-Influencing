# Title: Alpha diversity microbiome function
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

## Loading required package ----------------------------------------------------
library(tidyverse)
library(janitor)
library(knitr)
library(kableExtra)
library(vegan)
library(rstatix)
library(ggpubr)
library(robCompositions)
library(ape)
library(ggrepel)
library(patchwork)
library(ggprism)
library(ggnewscale)
## Setting color parameters ----------------------------------------------------
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02", "CI" = "#e7298a", 
               "SW" = "#b10c0c", "CDA" = "#ab37c9", "BI" = "#e96e0f",
               "TPC" = "#019a5b"
)
envPalette <- c("High" = "#D55E00","Low" = "#0072B2")
## Load raw data ---------------------------------------------------------------
Functional_list <- readRDS("~/Desktop/marine_sediment/results/functional_profile.rds")
## Functions -------------------------------------------------------------------
Alpha_diversity <- function(data,name){
  # data <- Functional_list$KO_list$KO_cpm
  # name <- "KO"
  
  taxon1 <- data %>%
    mutate(across(where(is.numeric), ~ as.integer(.))) %>%
    column_to_rownames("Feature") %>%
    t()
  
  richness_df <- estimateR(taxon1) %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::select(S.obs,S.chao1,S.ACE) %>%
    `colnames<-`(c("Richness","Chao1","ACE")) %>%
    rownames_to_column("SampleID")
  
  taxon <- data %>%
    column_to_rownames("Feature") %>%
    t()
  
  shannon_diversity <- as.data.frame(diversity(taxon)) %>%
    `colnames<-`("Shannon_diversity")
  simpson_diversity <- as.data.frame(diversity(taxon,index="simpson")) %>%
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
    pivot_longer(Shannon_diversity:ACE) %>%
    group_by(name) %>%
    wilcox_test(value ~ Gradient)
  
  # The palette with grey:
  cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02", "CI" = "#e7298a", 
                 "SW" = "#b10c0c", "CDA" = "#ab37c9", "BI" = "#e96e0f",
                 "TPC" = "#019a5b"
  )
  
  alpha_diversity_season <- alpha_diversity %>%
    pivot_longer(Shannon_diversity:ACE) %>%
    group_by(Site,Season,name) %>%
    summarise_at("value",mean) %>%
    ungroup() %>%
    group_by(name) %>%
    wilcox_test(value ~ Season,paired = T)
  
  alpha_diversity_all <- alpha_diversity %>%
    pivot_longer(Shannon_diversity:ACE) %>%
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
  
  # tool2 = alpha_diversity %>%
  #   wilcox_test(Shannon_diversity ~ Site) %>% 
  #   adjust_pvalue(method = 'fdr') %>%
  #   add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
  #                    symbols = c("****", "***", "**", "*", "ns")
  #   )
  gradient_p <- alpha_diversity_gradient %>% filter(name %in% "Shannon_diversity") %>%
    .$p %>% round(6)
  season_p <- alpha_diversity_season %>% filter(name %in% "Shannon_diversity") %>%
    .$p %>% round(2)
  gradient.label <- paste("Pollution P-value:",gradient_p)
  season.label <- bquote("Season P-value:"~.(season_p))
  nonmetric_label = list(gradient.label,season.label)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Shannon_diversity),0.999*max(alpha_diversity$Shannon_diversity))
  
  p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Shannon_diversity",outlier.shape=NA) +
    ylab(paste("Shannon diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Shannon_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    theme_prism(base_size = 16)
  #p1
  
  gradient_p <- alpha_diversity_gradient %>% filter(name %in% "Simpson_diversity") %>%
    .$p %>% round(6)
  season_p <- alpha_diversity_season %>% filter(name %in% "Simpson_diversity") %>%
    .$p %>% round(2)
  gradient.label <- paste("Pollution P-value:",gradient_p)
  season.label <- bquote("Season P-value:"~.(season_p))
  nonmetric_label = list(gradient.label,season.label)
  coord_x <- Inf
  coord_y <- c(0.9991,0.999)
  
  p2 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Simpson_diversity",outlier.shape=NA) +
    ylab(paste("Simpson diversity ","(",name,")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Simpson_diversity, shape=Season,color = Site), 
               position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    #annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
    theme_prism(base_size = 16)
  
  
  return_result <- list(Shannon=p1,
                        Simpson = p2,
                        all_p = p0,
                        sign_gradient = alpha_diversity_gradient,
                        sign_season = alpha_diversity_season,
                        alpha_diversity = alpha_diversity)
  return(return_result)
}
### Main part -----------------------------------------------------------------

COG_alphaD_list <- Alpha_diversity(Functional_list$COG_list$COG_cpm,"COG")
  

pdf("~/Desktop/marine_sediment/results/Function/Functional_alpha_diversity_COG.pdf",
    width=12, height= 6)
COG_alphaD_list$Shannon + COG_alphaD_list$Simpson +
  plot_layout(guides = 'collect')
dev.off()


return_list <- list(
  COG_alphaD_list = COG_alphaD_list
                    )
saveRDS(return_list, file="~/Desktop/marine_sediment/results/Function/functional_alpha_diversity.rds")



