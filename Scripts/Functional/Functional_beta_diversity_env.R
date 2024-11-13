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
library(patchwork)
## Setting color parameters ----------------------------------------------------
cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02", "CI" = "#e7298a", 
               "SW" = "#b10c0c", "CDA" = "#ab37c9", "BI" = "#e96e0f",
               "TPC" = "#019a5b")
envPalette <- c("High" = "#D55E00","Low" = "#0072B2")
## Load raw data ---------------------------------------------------------------
Prof_re <- readRDS("~/Desktop/marine_sediment/results/Functional_beta_diversity.rds")
## Significant species ---------------------------------------------------------
gradient_list <- readRDS("~/Desktop/marine_sediment/results/MEGAN6/Microbiome_DA_pollution_List.rds")
taxonomy_list <- list(
  phylum = gradient_list$phylum_list$taxon %>%
    as.data.frame() %>%
    dplyr::select(gradient_list$phylum_relist$sign_df$feature) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-1,names_to = "ID") %>%
    left_join(gradient_list$phylum_list$Feature_ID) %>%
    mutate(name = str_remove_all(name,"^.*?__")),
  genus = gradient_list$genus_list$taxon %>%
    as.data.frame() %>%
    dplyr::select(gradient_list$genus_relist$sign_df$feature) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-1,names_to = "ID") %>%
    left_join(gradient_list$genus_list$Feature_ID) %>%
    mutate(name = str_remove_all(name,"^.*?__")),
  species = gradient_list$species_list$taxon %>%
    as.data.frame() %>%
    dplyr::select(gradient_list$species_relist$sign_df$feature) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-1,names_to = "ID") %>%
    left_join(gradient_list$species_list$Feature_ID) %>%
    mutate(name = str_remove_all(name,"^.*?__"))
)
## Functions -------------------------------------------------------------------
func_envfit <- function(pca_tmp,microbiome_data){
  microbiome_data <- microbiome_data %>%
    dplyr::select(SampleID,name,value) %>%
    pivot_wider() %>%
    column_to_rownames("SampleID")
  fit <- envfit(pca_tmp,microbiome_data,perm=999)
  fit.df<-rbind(
    as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))
  )
  fit.df$name<-rownames(fit.df)
  fit.df$Rsquare <- c(fit$vectors$r)
  fit.df$pvals <- c(fit$vectors$pvals)
  fit.df$pvals.BH <- p.adjust(fit.df$pvals,method = "BH")
  fit.df.1 <- fit.df %>% 
    filter(pvals <= 0.05, pvals.BH <= 0.05)
  arrow_factor <- ordiArrowMul(fit)
  return_list <- list(
    df = fit.df.1,
    all_df = fit.df,
    arrow_factor = arrow_factor
  )
}
Plot_pcoa <- function(pcoa_list,adonis,name) {
  # pcoa_list <- Prof_re$LogEuclidean_pcoa[[1]]
  # adonis <- Prof_re$LogEuclidean_adonis[[1]]
  # name <- "Euclidean distance"
  
  Meta_reorder <- data.frame(ID = rownames(pcoa_list$points)) %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others)
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  title = sprintf("Principal Coordinates Analysis (%s)", name)
  eigen_percent <- round(((pcoa_list$eig)/sum(pcoa_list$eig))*100,2)
  pc1 <- paste("PCoA1 (",round(eigen_percent[1],digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(eigen_percent[2],digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=pcoa_list$points[,1],PC2=pcoa_list$points[,2],ID=rownames(pcoa_list$points)) %>%
    mutate(ID = str_replace_all(ID, "\\.", "-")) %>%
    left_join(Meta_reorder) %>%
    left_join(gradient)
  pca_tmp <- data.frame(PC1=pcoa_list$points[,1],PC2=pcoa_list$points[,2])
  phylum_fit <- func_envfit(pca_tmp,taxonomy_list$phylum)
  genus_fit <- func_envfit(pca_tmp,taxonomy_list$genus)
  species_fit <- func_envfit(pca_tmp,taxonomy_list$species)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis$adonis_gradient_R2)~"P = "~.(adonis$adonis_gradient_P))
  line_2 = bquote("PERMANOVA season:" ~ R^2 ~"= "~.(adonis$adonis_season_R2)~"P = "~.(adonis$adonis_season_P))

  nonmetric_label = list(line_1,line_2)
  coord_x <- Inf
  coord_y <- c(max(pcoa$PC2),0.9*max(pcoa$PC2))
  
  nonmetric_label = list(line_1,line_2)
  phylum_fit_df <- phylum_fit$df %>%
    filter(name %in% c("Proteobacteria","Planctomycetes",
                       "Nitrospinae","Chloroflexi",
                       "Candidatus Tectomicrobia","Bacteroidetes"
    ))
  multiply_num = 15
  if(name == "COG Pathway") {multiply_num = 2}
  p <- ggplot(pcoa, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    geom_segment(data=phylum_fit_df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* multiply_num,y=0,yend=PC2* phylum_fit$arrow_factor* multiply_num,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit_df,aes(x=PC1* phylum_fit$arrow_factor* multiply_num ,y=PC2* phylum_fit$arrow_factor* multiply_num,label=name),
                    size=5,colour="#EEAD0E",min.segment.length = 0)+
    labs(x=pc1,y=pc2) + 
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    new_scale_color() +
    labs(title = name)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,color=Gradient),
                 level=0.68,alpha=1,
                 type = "norm")+
    scale_color_manual(name = "Gradient",values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  return(p)
}
## Results ---------------------------------------------------------------------
Prof_plot <- Prof_re %>%
  mutate(
    plot1_LogEuclidean = pmap(list(LogEuclidean_pcoa,LogEuclidean_adonis,Category), 
                             function(pcoa_list,adonis,name) 
                               Plot_pcoa(pcoa_list,adonis,name)
    ),
    plot1_Bray = pmap(list(Bray_pcoa,Bray_adonis,Category), 
                     function(pcoa_list,adonis,name) 
                       Plot_pcoa(pcoa_list,adonis,name)
    )
  )


pdf("~/Desktop/marine_sediment/results/Function/Functional_beta_diversity_LogEuclidean_COG.pdf",
    width=12, height= 6)
Prof_plot$plot1_LogEuclidean[[1]] + Prof_plot$plot1_LogEuclidean[[3]] +
  plot_layout(guides = 'collect')
dev.off()
