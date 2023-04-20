library(tidyverse)
library(patchwork)
library(ggprism)
ancom2_gradient_list <- readRDS("~/Desktop/marine_sediment/RDS_files/ancom2_gradient_list.rds")
ancom2_season_mle4_list <- readRDS("~/Desktop/marine_sediment/RDS_files/ancom2_season_mle4_list.rds")
ancom2_season_list <- ancom2_season_mle4_list
data_path <- "~/Desktop/marine_sediment/data/"
gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                   Gradient = c("High","High","High","Medium","Medium","Low","Low"))
cbPalette <- c("#0C75D1", "#E6AB02", "#E7298A", "#B10C0C", "#AB37C9", "#E96E0F", "#019A5B")
#envPalette <- rev(c("#D95F02","#E6AB02","#1B9E77"))
envPalette <- c("#0072B2","#E69F00","#D55E00")
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
Plot_gradient <- function(taxon,feature_df){
  #taxon <- genus_otu
  #feature_df <- ancom2_gradient_list$Genus$re_1$df
  feature <- feature_df %>%
    filter(!W %in% "Inf") %>%
    filter(detected_0.7 %in% "TRUE") %>% 
    .$taxa_id
  RelativeAbundance_df <- taxon %>%
    mutate(tmp = SampleID) %>%
    separate(tmp,c("Site", "Season","tmp")) %>%
    dplyr::select(.,-tmp) %>%
    relocate(c("Site", "Season"),.after = SampleID) %>%
    as_tibble() %>%
    mutate_all(type.convert) %>%
    pivot_longer(4:(length(.)),names_to = "microbiome",values_to = "relativeAbundance") %>%
    filter(microbiome %in% feature) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient,levels = c("Low","Medium","High")))
  p <- ggplot(RelativeAbundance_df) + geom_boxplot(aes(x = Gradient, y = relativeAbundance,fill = Gradient),alpha = 0.2) +
    facet_wrap(~ microbiome,scales = "free") +
    geom_point(aes(x = Gradient, y = relativeAbundance,color = Site,shape = Season),size=3) +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = envPalette)+
    scale_shape_manual(values = c(17,16)) +
    theme_prism(base_size = 24) + 
    xlab("") + ylab("")+ 
    #theme_classic()+ 
    theme(legend.title = element_blank()) +
    theme(strip.background = element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(legend.text=element_text(size=24,face = "bold"))
  
  return(p)
}
Plot_season <- function(taxon,feature_df){
  # taxon <- genus_otu
  # feature_df <- ancom2_season_list$Genus$re_1$df
  feature <- feature_df %>%
    filter(!W %in% "Inf") %>%
    filter(detected_0.7 %in% "TRUE") %>% 
    .$taxa_id
  RelativeAbundance_df <- taxon %>%
    mutate(tmp = SampleID) %>%
    separate(tmp,c("Site", "Season","tmp")) %>%
    dplyr::select(.,-tmp) %>%
    relocate(c("Site", "Season"),.after = SampleID) %>%
    as_tibble() %>%
    mutate_all(type.convert) %>%
    pivot_longer(4:(length(.)),names_to = "microbiome",values_to = "relativeAbundance") %>%
    group_by(Site,Season,microbiome) %>%
    summarise_at("relativeAbundance", mean, na.rm = TRUE) %>%
    filter(microbiome %in% feature) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient,levels = c("Low","Medium","High"))) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter")))
  p <- ggplot(RelativeAbundance_df, aes(x=Season,y=relativeAbundance)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ microbiome,scales = "free") +
    geom_point(aes(color = Site),size=2, alpha=0.5) +
    geom_line(aes(group=Site,color=Site),linetype="11") +
    scale_colour_manual(values = cbPalette) +
    xlab("") +
    ylab("")+
    theme_prism()
    #theme_classic()
  return(p)
}
species_p <- Plot_gradient(species_otu,ancom2_gradient_list$Species$re_1$df)
species_p$facet$params$ncol <- 7
genus_p <- Plot_gradient(genus_otu,ancom2_gradient_list$Genus$re_1$df)
genus_p$facet$params$ncol <- 7
pdf("~/Desktop/marine_sediment/manuscript/Figure_v5/Tax_gradient_ordian.pdf",
    width = 32,height = 16)
genus_p/species_p +
  plot_layout(heights = c(1, 1),guides = 'collect')  &
  theme(legend.position='bottom')
dev.off()
genus_p <- Plot_gradient(genus_otu,ancom2_gradient_list$Genus$re_1$df)
genus_p$facet$params$ncol <- 7
order_p <- Plot_gradient(order_otu,ancom2_gradient_list$Order$re_1$df)
order_p$facet$params$ncol <- 4
phylum_p <- Plot_gradient(phylum_otu,ancom2_gradient_list$Phylum$re_1$df)
pdf("~/Desktop/marine_sediment/manuscript/Figure_v1/Tax_gradient_ordian.pdf",
    width = 14,height = 8)
(phylum_p + order_p + plot_layout(widths = c(2, 4))) / genus_p +
  plot_layout(heights = c(1, 3),guides = 'collect')  &
  theme(legend.position='bottom')
dev.off()
genus_p_season <- Plot_season(genus_otu,ancom2_season_list$Genus$re_1$df)
phylum_p_season <- Plot_season(phylum_otu,ancom2_season_list$Phylum$re_1$df)
Order_p_season <- Plot_season(order_otu,ancom2_season_list$Order$re_1$df)
genus_p_season$facet$params$ncol <- 7
pdf("~/Desktop/marine_sediment/manuscript/Figure_v1/Tax_season_ordian.pdf",width = 14,height = 5)
((phylum_p_season  + Order_p_season  + plot_spacer()  +plot_spacer()  +plot_spacer()  + plot_spacer() +  plot_layout(nrow = 1))/ genus_p_season) +
  plot_layout(heights = c(1, 1),guides = 'collect')  &
  theme(legend.position='bottom')
dev.off()              
