Functional_alpha_beta_diversity <- function()
{
### Loading required package ---------------------------------------------------
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
data_path <- "~/Desktop/marine_sediment/data/"
ancom2_gradient_list <- readRDS("~/Desktop/marine_sediment/RDS_files/ancom2_gradient_list.rds")
### Setting color parameters ---------------------------------------------------
cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
envPalette <- rev(c("#D55E00","#E69F00","#0072B2"))
### Load raw data --------------------------------------------------------------
# Assembly prok
COG_pathway <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_pathway_stat.flt.txt"
COG <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_cog.tsv"
COG_path_from_web <- read_tsv("~/Desktop/marine_sediment/results_new/functional_analyses/COG_pathway_from_web.tsv")
COG_pathway_df_raw <- read_tsv(COG_pathway) %>%
  filter(COG_pathway %in% COG_path_from_web$Pathway)
COG_df_raw <- read_tsv(COG)
COG_file <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_cog.flt.hellinger.tsv"
COG_path_file <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_pathway_stat.flt.hellinger.tsv"
### Read File ------------------------------------------------------------------
COG <- read_tsv(COG_file)
COG_path <- read_tsv(COG_path_file)

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

tidy_data <- function(data,sign)
{
  # data <- phylum_otu
  # sign <- ancom2_gradient_list$Phylum$re_1$df$taxa_id
  microbe <- data %>%
    dplyr::select(SampleID,sign) %>%
    mutate(tmp = SampleID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others,-Site) %>%
    column_to_rownames("SampleID")
    return(microbe) 
}

Phylum_list <- tidy_data(phylum_otu,ancom2_gradient_list$Phylum$re_1$df %>% 
                           filter(W != "Inf") %>%
                           filter(detected_0.7 %in% "TRUE") %>%
                           .$taxa_id)
Genus_list <- tidy_data(genus_otu,ancom2_gradient_list$Genus$re_1$df %>% 
                          filter(W != "Inf") %>%
                          filter(detected_0.7 %in% "TRUE") %>%
                          .$taxa_id)
Species_list <- tidy_data(species_otu,ancom2_gradient_list$Species$re_1$df %>% 
                            filter(W != "Inf") %>% 
                            filter(detected_0.7 %in% "TRUE") %>%
                            .$taxa_id)
taxonomy_list <- list(
  phylum = Phylum_list,
  genus = Genus_list,
  species = Species_list
    )
### Function Diversity ---------------------------------------------------------
# data <- COG_assembly_prok
# name <- "COG_assembly_prok"
Diversity_1 <- function(data1,data2,name,taxonomy_list)
{
  # data1 <- COG_pathway_df_raw
  # data2 <- COG_path
  # taxonomy_list <- taxonomy_list
  # name <- "COG_path_assembly_prok"
  ## Tidy data -----------------------------------------------------------------
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  taxon <- data1 %>%
    textshape::column_to_rownames(loc = 1) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_at(2:(length(.)),type.convert) %>%
    column_to_rownames() %>%
    as.data.frame() %>%
    replace(is.na(.), 0)
  
  # ## 50% prevalence filter in at leat one site ---------------------------------
  # Prev_df <- taxon %>%
  #   rownames_to_column() %>%
  #   mutate(sampleID = rowname) %>%
  #   separate(sampleID, into = c("Site", "Season","others"), sep = "_") %>%
  #   select(.,-others,-Season) %>%
  #   pivot_longer(2:(length(.)-1),names_to = "microbiome",values_to = "relativeAbundance") %>%
  #   group_by(microbiome,Site) %>%
  #   summarise_at("relativeAbundance", function(x) length(which(x!=0))) %>%
  #   ungroup() %>%
  #   group_by(microbiome) %>%
  #   summarise_at("relativeAbundance", function(x) length(which(x>=3))) %>%
  #   filter(relativeAbundance > 0)
  # taxon <- taxon %>%
  #   rownames_to_column() %>%
  #   pivot_longer(2:length(.),names_to = "microbiome",values_to = "relativeAbundance") %>%
  #   filter(microbiome %in% Prev_df$microbiome) %>%
  #   pivot_wider(names_from = microbiome,values_from = relativeAbundance) %>%
  #   column_to_rownames()
  
  shannon_diversity <- as.data.frame(diversity(taxon)) %>%
    `colnames<-`("Shannon_diversity")
  
  simpson_diversity <- as.data.frame(diversity(taxon,index="simpson")) %>%
    `colnames<-`("Simpson_diversity")
  
  alpha_diversity <- merge(shannon_diversity,simpson_diversity,by="row.names",
                           all.x = TRUE) %>%
    rename(sample = Row.names) %>%
    mutate(sampleID = sample) %>%
    separate(sampleID, into = c("Site", "Season","others"), sep = "_") %>%
    select(.,-others) %>%
    mutate(Season = str_replace_all(Season,c("S"="Summer","W"="Winter"))) %>%
    left_join(gradient)
  
  #### Figure 1 ----------------------------------------------------------------
  kruskal.p <- kruskal.test(Shannon_diversity~Gradient,data = alpha_diversity) %>% 
    tidy() %>% .$p.value %>% formatC(., format = "e", digits = 2)
  kruskal.p_tmp_list <- str_split(kruskal.p,"e")
  kruskal.p_tmp_list[[1]][2] = str_replace(kruskal.p_tmp_list[[1]][2],"^-0","-")
  kruskal.label <- bquote("Kruskal-Wallis H test P value: " ~.(kruskal.p_tmp_list[[1]][1])~"x"~10^~paste(.(kruskal.p_tmp_list[[1]][2])))
  #kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
  season_df <- alpha_diversity %>%
    group_by(Site,Season) %>%
    summarise_at("Shannon_diversity",mean) %>%
    ungroup()
  season_P <- wilcox.test(season_df %>% 
                            filter(Season %in% "Summer") %>%
                            .$Shannon_diversity,
                          season_df %>% 
                            filter(Season %in% "Winter") %>%
                            .$Shannon_diversity,
                          paired = T
  ) %>% tidy() %>% .$p.value %>% round(2)
  season.label <- bquote("Wilcoxon Signed-Rank Test P value:"~.(season_P))
  tool2 = alpha_diversity %>%
    wilcox_test(Shannon_diversity ~ Gradient) %>% 
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")
    )
  
  nonmetric_label = list(kruskal.label,season.label)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Shannon_diversity)  + 0.2 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)),
                                                                max(alpha_diversity$Shannon_diversity)  + 0.1 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)))
  
  alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
  p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Shannon_diversity",outlier.shape = NA
  ) + ylab(paste("Shannon diversity"," (","COG pathway",")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Shannon_diversity, shape=Season,color = Site), position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,size = 5) +
    # annotate("text", x = Inf, y =max(alpha_diversity$Shannon_diversity) + 0.2 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)),
    #          label= kruskal.label,hjust = 1, vjust = 1,size = 6) + 
    stat_pvalue_manual(tool2, label = "p.adj.signif", 
                       hide.ns = TRUE,
                       size = 6,
                       y.position = seq(max(alpha_diversity$Shannon_diversity),
                                        max(alpha_diversity$Shannon_diversity)+0.05*length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]),
                                        length.out=length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]))
    ) +
    theme_prism(base_size = 12)
  #theme(legend.position="none")
  
  #### Figure 2 ----------------------------------------------------------------
  kruskal.p <- kruskal.test(Simpson_diversity~Gradient,data = alpha_diversity) %>% 
    tidy() %>% .$p.value %>% formatC(., format = "e", digits = 2)
  kruskal.p_tmp_list <- str_split(kruskal.p,"e")
  kruskal.p_tmp_list[[1]][2] = str_replace(kruskal.p_tmp_list[[1]][2],"^-0","-")
  kruskal.label <- bquote("Kruskal-Wallis H test P value: " ~.(kruskal.p_tmp_list[[1]][1])~"x"~10^~paste(.(kruskal.p_tmp_list[[1]][2])))
  #kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
  season_df <- alpha_diversity %>%
    group_by(Site,Season) %>%
    summarise_at("Simpson_diversity",mean) %>%
    ungroup()
  season_P <- wilcox.test(season_df %>% 
                            filter(Season %in% "Summer") %>%
                            .$Simpson_diversity,
                          season_df %>% 
                            filter(Season %in% "Winter") %>%
                            .$Simpson_diversity,
                          paired = T
  ) %>% tidy() %>% .$p.value %>% round(2)
  season.label <- bquote("Wilcoxon Signed-Rank Test P value:"~.(season_P))
  tool3 = alpha_diversity %>%
    wilcox_test(Simpson_diversity ~ Gradient) %>% 
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")
    )
  
  nonmetric_label = list(kruskal.label,season.label)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Simpson_diversity)  + 0.2 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)),
               max(alpha_diversity$Simpson_diversity)  + 0.1 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)))
  
  #alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  #alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
  p2 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Simpson_diversity",outlier.shape = NA
  ) + ylab(paste("Simpson diversity"," (","COG pathway",")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Simpson_diversity, shape=Season,color = Site), position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,size = 5) +
    # annotate("text", x = Inf, y =max(alpha_diversity$Simpson_diversity) + 0.2 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)),
    #          label= kruskal.label,hjust = 1, vjust = 1,size = 6) + 
    stat_pvalue_manual(tool3, label = "p.adj.signif", 
                       hide.ns = TRUE,
                       size = 6,
                       y.position = seq(max(alpha_diversity$Simpson_diversity),
                                        max(alpha_diversity$Simpson_diversity)+0.000008*length(tool3$p.adj.signif[tool3$p.adj.signif != "ns"]),
                                        length.out=length(tool3$p.adj.signif[tool3$p.adj.signif != "ns"]))
    ) +
    theme_prism(base_size = 12)
  
  
  ### beta-diversity ----------------------------------------------------------------
  taxon <- data2 %>%
    textshape::column_to_rownames(loc = 1) %>%
    as.matrix() %>%
    #t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_at(2:(length(.)),type.convert) %>%
    column_to_rownames() %>%
    as.data.frame() %>%
    replace(is.na(.), 0)
  ####bray_curtis
  taxon <- as.matrix(taxon)
  bray_curtis <- vegdist(taxon, method = "bray")
  bray_curtis <- as.matrix(bray_curtis)
  ####Artchison
  min_tmp <- min(taxon[taxon>0])
  taxon[taxon == 0] <- min_tmp
  Artchison_dist <- as.matrix(aDist(taxon))
  ####Euclidean distance
  Euclidean_dist <- vegdist(taxon, method = "euclidean")
  Euclidean_dist <- as.matrix(Euclidean_dist)
  ###PERMANOVA
  adonis_Euclidean_dist = adonis2(Euclidean_dist~Gradient, alpha_diversity, permutations = 999)
  adonis_Euclidean_R2 = round(adonis_Euclidean_dist$R2[1],2)
  adonis_Euclidean_P = round(adonis_Euclidean_dist$`Pr(>F)`[1],4)
  
  adonis_Euclidean_dist_1 = adonis2(Euclidean_dist~Season, alpha_diversity, permutations = 999)
  adonis_Euclidean_R2_1 = round(adonis_Euclidean_dist_1$R2[1],2)
  adonis_Euclidean_P_1 = round(adonis_Euclidean_dist_1$`Pr(>F)`[1],4)
  
  adonis_bray_curtis_dist = adonis2(bray_curtis~Gradient, alpha_diversity, permutations = 999)
  adonis_bray_curtis_R2 = round(adonis_bray_curtis_dist$R2[1],2)
  adonis_bray_curtis_P = round(adonis_bray_curtis_dist$`Pr(>F)`[1],4)
  
  adonis_bray_curtis_dist_1 = adonis2(bray_curtis~Season, alpha_diversity, permutations = 999)
  adonis_bray_curtis_R2_1 = round(adonis_bray_curtis_dist_1$R2[1],2)
  adonis_bray_curtis_P_1 = round(adonis_bray_curtis_dist_1$`Pr(>F)`[1],4)
  
  adonis_Artchison_dist = adonis2(Artchison_dist~Gradient, alpha_diversity, permutations = 999)
  adonis_Artchison_R2 = round(adonis_Artchison_dist$R2[1],2)
  adonis_Artchison_P = round(adonis_Artchison_dist$`Pr(>F)`[1],4)
  
  adonis_Artchison_dist_1 = adonis2(Artchison_dist~Season, alpha_diversity, permutations = 999)
  adonis_Artchison_R2_1 = round(adonis_Artchison_dist_1$R2[1],2)
  adonis_Artchison_P_1 = round(adonis_Artchison_dist_1$`Pr(>F)`[1],4)
  ###betadisper
  beta_Euclidean_dist <- betadisper(as.dist(Euclidean_dist), alpha_diversity$Gradient) %>% permutest()
  beta_Euclidean_P <- round(as.data.frame(beta_Euclidean_dist$tab)[1,]$`Pr(>F)`,4)
  
  beta_bray_curtis_dist <- betadisper(as.dist(bray_curtis), alpha_diversity$Gradient) %>% permutest()
  beta_bray_curtis_P <- round(as.data.frame(beta_bray_curtis_dist$tab)[1,]$`Pr(>F)`,4)
  
  beta_Artchison_dist <- betadisper(as.dist(Artchison_dist), alpha_diversity$Gradient) %>% permutest()
  beta_Artchison_P <- round(as.data.frame(beta_Artchison_dist$tab)[1,]$`Pr(>F)`,4)
  
  ### microbiome data
  phylum <- taxonomy_list$phylum
  genus <- taxonomy_list$genus
  species <- taxonomy_list$species
  func_envfit <- function(pca_tmp,microbiome_data){
    fit <- envfit(pca_tmp,microbiome_data,perm=999)
    fit.df<-rbind(
      as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r)),
      as.data.frame(fit$factors$centroids*sqrt(fit$factors$r))
    )
    fit.df$name<-rownames(fit.df)
    fit.df$Rsquare <- c(fit$vectors$r,
                            c("SeasonSummer" = fit$factors$r[[1]],
                              "SeasonWinter" = fit$factors$r[[1]]))
    fit.df$pvals <- c(fit$vectors$pvals,
                      c("SeasonSummer" = fit$factors$pvals[[1]],
                        "SeasonWinter" = fit$factors$pvals[[1]]))
    fit.df$pvals.BH <- p.adjust(fit.df$pvals,method = "BH")
    fit.df.1 <- fit.df %>% 
      filter(pvals < 0.05, pvals.BH < 0.1)
    fit.df.2 <- fit.df %>% 
      filter(pvals.BH >= 0.1)
    arrow_factor <- ordiArrowMul(fit)
    return_list <- list(
      df = fit.df.1,
      tmp = fit.df,
      arrow_factor = arrow_factor
    )
  }
  ####PcoA -----------------------------------------------------------------------
  #####bray_curtis
  PCOA <- pcoa(bray_curtis, correction="none", rn=NULL)
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis (%s) (%s)",name,"bray curtis distance")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  pcoa_tmp <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2])
  phylum_fit <- func_envfit(pcoa_tmp,phylum)
  genus_fit <- func_envfit(pcoa_tmp,genus)
  species_fit <- func_envfit(pcoa_tmp,species)
  line_1 = bquote("PERMANOVA gradient:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2)~"P = "~.(adonis_bray_curtis_P))
  line_2 = bquote("PERMDISP2 gradient:"~"P = "~.(beta_bray_curtis_P))
  line_3 = bquote("PERMANOVA season:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2_1)~"P = "~.(adonis_bray_curtis_P_1))
  nonmetric_label = list(line_1,line_2,line_3)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(pcoa$PC2),0.9*max(pcoa$PC2),0.8*max(pcoa$PC2))
  
  p3 <- ggplot(pcoa, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) + 
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    new_scale_color() +
    labs(title = "COG pathway")+
    geom_segment(data=phylum_fit$df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* 0.01,y=0,yend=PC2* phylum_fit$arrow_factor* 0.01,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit$df,aes(x=PC1* phylum_fit$arrow_factor* 0.01 ,y=PC2* phylum_fit$arrow_factor* 0.01,label=name),
                    size=3,colour="#EEAD0E",min.segment.length = 0.01)+
    geom_segment(data=genus_fit$df,aes(x=0,xend=PC1* genus_fit$arrow_factor* 0.01,y=0,yend=PC2* genus_fit$arrow_factor* 0.01,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    geom_text_repel(data=genus_fit$df,aes(x=PC1* genus_fit$arrow_factor * 0.01,y=PC2* genus_fit$arrow_factor* 0.01,label=name),
                    size=3,colour="#6495ED",min.segment.length = 0.01)+
    geom_segment(data=species_fit$df,aes(x=0,xend=PC1* species_fit$arrow_factor * 0.01,y=0,yend=PC2* species_fit$arrow_factor * 0.01),
                 arrow = arrow(length = unit(0.2, "cm")),colour="grey",alpha = 0.5) +
    geom_text_repel(data=species_fit$df,aes(x=PC1* species_fit$arrow_factor * 0.01,y=PC2* species_fit$arrow_factor * 0.01,label=name),
                    size=3,color = "grey",min.segment.length = 0.01) +
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,color=Gradient),
                 level=0.68,alpha=1,
                 type = "norm")+
    scale_color_manual(name = "Gradient",values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #####Artchison ------
  PCOA <- pcoa(Artchison_dist, correction="none", rn=NULL)
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis (%s) (%s)",name,"Artchison distance")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  pcoa_tmp <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2])
  phylum_fit <- func_envfit(pcoa_tmp,phylum)
  genus_fit <- func_envfit(pcoa_tmp,genus)
  species_fit <- func_envfit(pcoa_tmp,species)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Artchison_R2)~"P = "~.(adonis_Artchison_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Artchison_P))
  nonmetric_label = list(line_1,line_2)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(pcoa$PC2),0.80*max(pcoa$PC2))
  
  p4 <- ggplot(pcoa, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) + 
    geom_segment(data=phylum_fit$df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* 1,y=0,yend=PC2* phylum_fit$arrow_factor* 1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit$df,aes(x=PC1* phylum_fit$arrow_factor* 1 ,y=PC2* phylum_fit$arrow_factor* 1,label=name),
                    size=3,colour="#EEAD0E",min.segment.length = 0)+
    geom_segment(data=genus_fit$df,aes(x=0,xend=PC1* genus_fit$arrow_factor* 1,y=0,yend=PC2* genus_fit$arrow_factor* 1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    geom_text_repel(data=genus_fit$df,aes(x=PC1* genus_fit$arrow_factor * 1,y=PC2* genus_fit$arrow_factor* 1,label=name),
                    size=3,colour="#6495ED",min.segment.length = 0)+
    geom_segment(data=species_fit$df,aes(x=0,xend=PC1* species_fit$arrow_factor * 1,y=0,yend=PC2* species_fit$arrow_factor * 1),
                 arrow = arrow(length = unit(0.2, "cm")),colour="grey",alpha = 0.5) +
    geom_text_repel(data=species_fit$df,aes(x=PC1* species_fit$arrow_factor * 1,y=PC2* species_fit$arrow_factor * 1,label=name),
                    size=3,color = "grey",min.segment.length = 0) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  ####NMDS -----------------------------------------------------------------------
  #####bray_curtis
  sol <- metaMDS(bray_curtis)
  stress  <- as.numeric(sprintf("%.2f",sol$stress))
  title = sprintf("Non-metric Multidimensional Scaling (%s) (%s)",name,"bray curtis distance")
  pc1 <- "NMDS 1"
  pc2 <- "NMDS 2"
  nmds <- data.frame(PC1=sol$points[,1],PC2=sol$points[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2)~"P = "~.(adonis_bray_curtis_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_bray_curtis_P))
  line_3 = bquote("Stress:"~.(stress))
  nonmetric_label = list(line_1,line_2,line_3)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(nmds$PC2),0.90*max(nmds$PC2),0.80*max(nmds$PC2))
  
  p5 <- ggplot(nmds, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) +
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #####Artchison
  sol <- metaMDS(Artchison_dist)
  stress  <- as.numeric(sprintf("%.2f",sol$stress))
  title = sprintf("Non-metric Multidimensional Scaling (%s) (%s)",name,"Artchison distance")
  pc1 <- "NMDS 1"
  pc2 <- "NMDS 2"
  nmds <- data.frame(PC1=sol$points[,1],PC2=sol$points[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Artchison_R2)~"P = "~.(adonis_Artchison_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Artchison_P))
  line_3 = bquote("Stress:"~.(stress))
  nonmetric_label = list(line_1,line_2,line_3)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(nmds$PC2),0.85*max(nmds$PC2),0.70*max(nmds$PC2))
  
  p6 <- ggplot(nmds, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) +
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #final_plot <- (p1 + p2)/(p3 + p4) + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A')
  return_list <- list(P_Shannon = p1,
                      P_Simpson = p2,
                      #P_pca = pca_plot,
                      P_PCoA_bray = p3,
                      P_PCoA_Artchison = p4,
                      P_NMDS_bray = p5,
                      P_NMDS_Artchison = p6
  )
  return(return_list)
} 
Diversity_2 <- function(data1,data2,name,taxonomy_list)
{
  # data1 <- COG_df_raw
  # data2 <- COG
  # taxonomy_list <- taxonomy_list
  # name <- "COG_assembly_prok"
  ## Tidy data -----------------------------------------------------------------
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  taxon <- data1 %>%
    textshape::column_to_rownames(loc = 1) %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_at(2:(length(.)),type.convert) %>%
    column_to_rownames() %>%
    as.data.frame() %>%
    replace(is.na(.), 0)
  
  # ## 50% prevalence filter in at leat one site ---------------------------------
  # Prev_df <- taxon %>%
  #   rownames_to_column() %>%
  #   mutate(sampleID = rowname) %>%
  #   separate(sampleID, into = c("Site", "Season","others"), sep = "_") %>%
  #   select(.,-others,-Season) %>%
  #   pivot_longer(2:(length(.)-1),names_to = "microbiome",values_to = "relativeAbundance") %>%
  #   group_by(microbiome,Site) %>%
  #   summarise_at("relativeAbundance", function(x) length(which(x!=0))) %>%
  #   ungroup() %>%
  #   group_by(microbiome) %>%
  #   summarise_at("relativeAbundance", function(x) length(which(x>=3))) %>%
  #   filter(relativeAbundance > 0)
  # taxon <- taxon %>%
  #   rownames_to_column() %>%
  #   pivot_longer(2:length(.),names_to = "microbiome",values_to = "relativeAbundance") %>%
  #   filter(microbiome %in% Prev_df$microbiome) %>%
  #   pivot_wider(names_from = microbiome,values_from = relativeAbundance) %>%
  #   column_to_rownames()
  
  shannon_diversity <- as.data.frame(diversity(taxon)) %>%
    `colnames<-`("Shannon_diversity")
  
  simpson_diversity <- as.data.frame(diversity(taxon,index="simpson")) %>%
    `colnames<-`("Simpson_diversity")
  
  alpha_diversity <- merge(shannon_diversity,simpson_diversity,by="row.names",
                           all.x = TRUE) %>%
    rename(sample = Row.names) %>%
    mutate(sampleID = sample) %>%
    separate(sampleID, into = c("Site", "Season","others"), sep = "_") %>%
    select(.,-others) %>%
    mutate(Season = str_replace_all(Season,c("S"="Summer","W"="Winter"))) %>%
    left_join(gradient)
  
  #### Figure 1 ----------------------------------------------------------------
  kruskal.p <- kruskal.test(Shannon_diversity~Gradient,data = alpha_diversity) %>%
    tidy() %>% .$p.value  %>% formatC(., format = "e", digits = 2)
  kruskal.p_tmp_list <- str_split(kruskal.p,"e")
  kruskal.p_tmp_list[[1]][2] = str_replace(kruskal.p_tmp_list[[1]][2],"^-0","-")
  kruskal.label <- bquote("Kruskal-Wallis H test P value: " ~.(kruskal.p_tmp_list[[1]][1])~"x"~10^~paste(.(kruskal.p_tmp_list[[1]][2])))
  #kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
  season_df <- alpha_diversity %>%
    group_by(Site,Season) %>%
    summarise_at("Shannon_diversity",mean) %>%
    ungroup()
  season_P <- wilcox.test(season_df %>% 
                            filter(Season %in% "Summer") %>%
                            .$Shannon_diversity,
                          season_df %>% 
                            filter(Season %in% "Winter") %>%
                            .$Shannon_diversity,
                          paired = T
  ) %>% tidy() %>% .$p.value %>% round(2)
  #season.label <- bquote("Wilcoxon Signed-Rank Test P value:"~.(season_P))
  tool2 = alpha_diversity %>%
    wilcox_test(Shannon_diversity ~ Gradient) %>% 
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")
    )
  
  #nonmetric_label = list(kruskal.label,season.label)
  nonmetric_label = list(kruskal.label)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Shannon_diversity)  + 0.2 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)))
  
  # coord_y <- c(max(alpha_diversity$Shannon_diversity)  + 0.2 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)),
  #              max(alpha_diversity$Shannon_diversity)  + 0.1 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)))
  # 
  alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
  p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Shannon_diversity",outlier.shape = NA
  ) + ylab(paste("Shannon diversity"," (","COG term",")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Shannon_diversity, shape=Season,color = Site), position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,size = 5,
             size = 16/.pt,fontface = "bold") +
    # annotate("text", x = Inf, y =max(alpha_diversity$Shannon_diversity) + 0.2 * (max(alpha_diversity$Shannon_diversity) - min(alpha_diversity$Shannon_diversity)),
    #          label= kruskal.label,hjust = 1, vjust = 1,size = 6) + 
    stat_pvalue_manual(tool2, label = "p.adj.signif", 
                       size = 6,
                       hide.ns = TRUE,
                       y.position = seq(max(alpha_diversity$Shannon_diversity),
                                        max(alpha_diversity$Shannon_diversity)+0.05*length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]),
                                        length.out=length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]))
    ) +
    theme_prism(base_size = 16)
  #theme(legend.position="none")
  
  #### Figure 2 ----------------------------------------------------------------
  kruskal.p <- kruskal.test(Simpson_diversity~Gradient,data = alpha_diversity) %>% 
    tidy() %>% .$p.value %>% formatC(., format = "e", digits = 2)
  kruskal.p_tmp_list <- str_split(kruskal.p,"e")
  kruskal.p_tmp_list[[1]][2] = str_replace(kruskal.p_tmp_list[[1]][2],"^-0","-")
  kruskal.label <- bquote("Kruskal-Wallis H test P value: " ~.(kruskal.p_tmp_list[[1]][1])~"x"~10^~paste(.(kruskal.p_tmp_list[[1]][2])))
  #kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
  season_df <- alpha_diversity %>%
    group_by(Site,Season) %>%
    summarise_at("Simpson_diversity",mean) %>%
    ungroup()
  season_P <- wilcox.test(season_df %>% 
                            filter(Season %in% "Summer") %>%
                            .$Simpson_diversity,
                          season_df %>% 
                            filter(Season %in% "Winter") %>%
                            .$Simpson_diversity,
                          paired = T
  ) %>% tidy() %>% .$p.value %>% round(2)
  #season.label <- bquote("Wilcoxon Signed-Rank Test P value:"~.(season_P))
  tool3 = alpha_diversity %>%
    wilcox_test(Simpson_diversity ~ Gradient) %>% 
    adjust_pvalue(method = 'fdr') %>%
    add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", "ns")
    )
  
  #nonmetric_label = list(kruskal.label,season.label)
  nonmetric_label = list(kruskal.label)
  coord_x <- Inf
  coord_y <- c(max(alpha_diversity$Simpson_diversity)  + 0.3 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)))
  
  # coord_y <- c(max(alpha_diversity$Simpson_diversity)  + 0.3 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)),
  #              max(alpha_diversity$Simpson_diversity)  + 0.2 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)))
  # 
  #alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  #alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
  p2 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Simpson_diversity",outlier.shape = NA
  ) + ylab(paste("Simpson diversity"," (","COG term",")",sep = "")) +
    geom_point(data = alpha_diversity, aes(x = Gradient, y = Simpson_diversity, shape=Season,color = Site), position = position_jitter(width = 0.2))+
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) + xlab("")+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,size = 5,
             size = 16/.pt,fontface = "bold") +
    # annotate("text", x = Inf, y =max(alpha_diversity$Simpson_diversity) + 0.2 * (max(alpha_diversity$Simpson_diversity) - min(alpha_diversity$Simpson_diversity)),
    #          label= kruskal.label,hjust = 1, vjust = 1,size = 6) + 
    stat_pvalue_manual(tool3, label = "p.adj.signif", 
                       hide.ns = TRUE,
                       size = 6,
                       y.position = seq(max(alpha_diversity$Simpson_diversity),
                                        max(alpha_diversity$Simpson_diversity)+0.000008*length(tool3$p.adj.signif[tool3$p.adj.signif != "ns"]),
                                        length.out=length(tool3$p.adj.signif[tool3$p.adj.signif != "ns"]))
    ) +
    theme_prism(base_size = 16)
  
  ### beta-diversity ----------------------------------------------------------------
  taxon <- data2 %>%
    textshape::column_to_rownames(loc = 1) %>%
    as.matrix() %>%
    #t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate_at(2:(length(.)),type.convert) %>%
    column_to_rownames() %>%
    as.data.frame() %>%
    replace(is.na(.), 0)
  ####bray_curtis
  taxon <- as.matrix(taxon)
  bray_curtis <- vegdist(taxon, method = "bray")
  bray_curtis <- as.matrix(bray_curtis)
  ####Artchison
  min_tmp <- min(taxon[taxon>0])
  taxon[taxon == 0] <- min_tmp
  Artchison_dist <- as.matrix(aDist(taxon))
  ####Euclidean distance
  Euclidean_dist <- vegdist(taxon, method = "euclidean")
  Euclidean_dist <- as.matrix(Euclidean_dist)
  ###PERMANOVA
  adonis_Euclidean_dist = adonis2(Euclidean_dist~Gradient, alpha_diversity, permutations = 999)
  adonis_Euclidean_R2 = round(adonis_Euclidean_dist$R2[1],2)
  adonis_Euclidean_P = round(adonis_Euclidean_dist$`Pr(>F)`[1],4)
  
  adonis_Euclidean_dist_1 = adonis2(Euclidean_dist~Season, alpha_diversity, permutations = 999)
  adonis_Euclidean_R2_1 = round(adonis_Euclidean_dist_1$R2[1],2)
  adonis_Euclidean_P_1 = round(adonis_Euclidean_dist_1$`Pr(>F)`[1],4)
  
  adonis_bray_curtis_dist = adonis2(bray_curtis~Gradient, alpha_diversity, permutations = 999)
  adonis_bray_curtis_R2 = round(adonis_bray_curtis_dist$R2[1],2)
  adonis_bray_curtis_P = round(adonis_bray_curtis_dist$`Pr(>F)`[1],4)
  
  adonis_bray_curtis_dist_1 = adonis2(bray_curtis~Season, alpha_diversity, permutations = 999)
  adonis_bray_curtis_R2_1 = round(adonis_bray_curtis_dist_1$R2[1],2)
  adonis_bray_curtis_P_1 = round(adonis_bray_curtis_dist_1$`Pr(>F)`[1],4)
  
  adonis_Artchison_dist = adonis2(Artchison_dist~Gradient, alpha_diversity, permutations = 999)
  adonis_Artchison_R2 = round(adonis_Artchison_dist$R2[1],2)
  adonis_Artchison_P = round(adonis_Artchison_dist$`Pr(>F)`[1],4)
  
  adonis_Artchison_dist_1 = adonis2(Artchison_dist~Season, alpha_diversity, permutations = 999)
  adonis_Artchison_R2_1 = round(adonis_Artchison_dist_1$R2[1],2)
  adonis_Artchison_P_1 = round(adonis_Artchison_dist_1$`Pr(>F)`[1],4)
  ###betadisper
  beta_Euclidean_dist <- betadisper(as.dist(Euclidean_dist), alpha_diversity$Gradient) %>% permutest()
  beta_Euclidean_P <- round(as.data.frame(beta_Euclidean_dist$tab)[1,]$`Pr(>F)`,4)
  
  beta_bray_curtis_dist <- betadisper(as.dist(bray_curtis), alpha_diversity$Gradient) %>% permutest()
  beta_bray_curtis_P <- round(as.data.frame(beta_bray_curtis_dist$tab)[1,]$`Pr(>F)`,4)
  
  beta_Artchison_dist <- betadisper(as.dist(Artchison_dist), alpha_diversity$Gradient) %>% permutest()
  beta_Artchison_P <- round(as.data.frame(beta_Artchison_dist$tab)[1,]$`Pr(>F)`,4)
  
  ### microbiome data
  phylum <- taxonomy_list$phylum
  genus <- taxonomy_list$genus
  species <- taxonomy_list$species
  func_envfit <- function(pca_tmp,microbiome_data){
    fit <- envfit(pca_tmp,microbiome_data,perm=999)
    fit.df<-rbind(
      as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r)),
      as.data.frame(fit$factors$centroids*sqrt(fit$factors$r))
    )
    fit.df$name<-rownames(fit.df)
    fit.df$Rsquare <- c(fit$vectors$r,
                        c("SeasonSummer" = fit$factors$r[[1]],
                          "SeasonWinter" = fit$factors$r[[1]]))
    fit.df$pvals <- c(fit$vectors$pvals,
                      c("SeasonSummer" = fit$factors$pvals[[1]],
                        "SeasonWinter" = fit$factors$pvals[[1]]))
    fit.df$pvals.BH <- p.adjust(fit.df$pvals,method = "BH")
    fit.df.1 <- fit.df %>% 
      filter(pvals < 0.05, pvals.BH < 0.1)
    fit.df.2 <- fit.df %>% 
      filter(pvals.BH >= 0.1)
    arrow_factor <- ordiArrowMul(fit)
    return_list <- list(
      df = fit.df.1,
      tmp = fit.df,
      arrow_factor = arrow_factor
    )
  }
  ## PCA -----------------------------------------------------------------------
  pca.out <- prcomp(taxon)
  pr.var<-pca.out$sdev^2
  pve<-pr.var/sum(pr.var)
  title = sprintf("Principle Component Analysis (%s)",name)
  pc1 <- paste("PC1 (",round(pve[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PC2 (",round(pve[2]*100,digits=2),"%)",sep="")
  pca <- data.frame(PC1=pca.out$x[,1],PC2=pca.out$x[,2],ind=alpha_diversity$sample,
                    Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  pca_tmp <- data.frame(PC1=pca.out$x[,1],PC2=pca.out$x[,2])
  phylum_fit <- func_envfit(pca_tmp,phylum)
  genus_fit <- func_envfit(pca_tmp,genus)
  species_fit <- func_envfit(pca_tmp,species)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Euclidean_R2)~"P = "~.(adonis_Euclidean_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Euclidean_P))
  #line_3 = bquote("ANOSIM:"~"R"~"= "~.(anosim_Euclidean_R)~"P = "~.(anosim_Euclidean_P))
  nonmetric_label = list(line_1,line_2)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(pca$PC2),0.90*max(pca$PC2))
  pca_plot <- ggplot(pca, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) + 
    geom_segment(data=phylum_fit$df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* 0.1,y=0,yend=PC2* phylum_fit$arrow_factor* 0.1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit$df,aes(x=PC1* phylum_fit$arrow_factor* 0.1 ,y=PC2* phylum_fit$arrow_factor* 0.1,label=name),
                    size=5,colour="#EEAD0E",min.segment.length = 0)+
    geom_segment(data=genus_fit$df,aes(x=0,xend=PC1* genus_fit$arrow_factor* 0.1,y=0,yend=PC2* genus_fit$arrow_factor* 0.1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    geom_text_repel(data=genus_fit$df,aes(x=PC1* genus_fit$arrow_factor * 0.1,y=PC2* genus_fit$arrow_factor* 0.1,label=name),
                    size=5,colour="#6495ED",min.segment.length = 0)+
    geom_segment(data=species_fit$df,aes(x=0,xend=PC1* species_fit$arrow_factor * 0.1,y=0,yend=PC2* species_fit$arrow_factor * 0.1),
                 arrow = arrow(length = unit(0.2, "cm")),colour="black",alpha = 0.5) +
    geom_text_repel(data=species_fit$df,aes(x=PC1* species_fit$arrow_factor * 0.1,y=PC2* species_fit$arrow_factor * 0.1,label=name),
                    size=5,color = "black",min.segment.length = 0) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  ####PcoA -----------------------------------------------------------------------
  #####bray_curtis
  PCOA <- pcoa(bray_curtis, correction="none", rn=NULL)
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis (%s) (%s)",name,"bray curtis distance")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  pcoa_tmp <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2])
  phylum_fit <- func_envfit(pcoa_tmp,phylum)
  genus_fit <- func_envfit(pcoa_tmp,genus)
  species_fit <- func_envfit(pcoa_tmp,species)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2)~"P = "~.(adonis_bray_curtis_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_bray_curtis_P))
  #line_3 = bquote("PERMANOVA season:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2_1)~"P = "~.(adonis_bray_curtis_P_1))
  #nonmetric_label = list(line_1,line_2,line_3)
  nonmetric_label = list(line_1,line_2)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(1.1*max(pcoa$PC2),0.95*max(pcoa$PC2))
  #coord_y <- c(max(pcoa$PC2),0.9*max(pcoa$PC2),0.85*max(pcoa$PC2))
  
  p3 <- ggplot(pcoa, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) + 
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    new_scale_color() +
    labs(title = "COG term")+
    geom_segment(data=phylum_fit$df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* 0.03,y=0,yend=PC2* phylum_fit$arrow_factor* 0.03,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit$df,aes(x=PC1* phylum_fit$arrow_factor* 0.03 ,y=PC2* phylum_fit$arrow_factor* 0.03,label=name),
                    size=5,colour="#EEAD0E",min.segment.length = 0.03)+
    geom_segment(data=genus_fit$df,aes(x=0,xend=PC1* genus_fit$arrow_factor* 0.03,y=0,yend=PC2* genus_fit$arrow_factor* 0.03,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    geom_text_repel(data=genus_fit$df,aes(x=PC1* genus_fit$arrow_factor * 0.03,y=PC2* genus_fit$arrow_factor* 0.03,label=name),
                    size=5,colour="#6495ED",min.segment.length = 0.03)+
    geom_segment(data=species_fit$df,aes(x=0,xend=PC1* species_fit$arrow_factor * 0.03,y=0,yend=PC2* species_fit$arrow_factor * 0.03),
                 arrow = arrow(length = unit(0.2, "cm")),colour="black",alpha = 0.5) +
    geom_text_repel(data=species_fit$df,aes(x=PC1* species_fit$arrow_factor * 0.03,y=PC2* species_fit$arrow_factor * 0.03,label=name),
                    size=5,color = "black",min.segment.length = 0.03) +
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt,fontface = "bold") +
    theme_prism(base_size = 16,border = F)+
    stat_ellipse(aes(group = Gradient,color=Gradient),
                 level=0.68,alpha=1,
                 type = "norm")+
    scale_color_manual(name = "Gradient",values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #####Artchison
  PCOA <- pcoa(Artchison_dist, correction="none", rn=NULL)
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis (%s) (%s)",name,"Artchison distance")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  pcoa_tmp <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2])
  phylum_fit <- func_envfit(pcoa_tmp,phylum)
  genus_fit <- func_envfit(pcoa_tmp,genus)
  species_fit <- func_envfit(pcoa_tmp,species)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Artchison_R2)~"P = "~.(adonis_Artchison_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Artchison_P))
  nonmetric_label = list(line_1,line_2)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(pcoa$PC2),0.80*max(pcoa$PC2))
  
  p4 <- ggplot(pcoa, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) + 
    geom_segment(data=phylum_fit$df,aes(x=0,xend=PC1* phylum_fit$arrow_factor* 1,y=0,yend=PC2* phylum_fit$arrow_factor* 1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=phylum_fit$df,aes(x=PC1* phylum_fit$arrow_factor* 1 ,y=PC2* phylum_fit$arrow_factor* 1,label=name),
                    size=3,colour="#EEAD0E",min.segment.length = 0)+
    geom_segment(data=genus_fit$df,aes(x=0,xend=PC1* genus_fit$arrow_factor* 1,y=0,yend=PC2* genus_fit$arrow_factor* 1,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    geom_text_repel(data=genus_fit$df,aes(x=PC1* genus_fit$arrow_factor * 1,y=PC2* genus_fit$arrow_factor* 1,label=name),
                    size=3,colour="#6495ED",min.segment.length = 0)+
    geom_segment(data=species_fit$df,aes(x=0,xend=PC1* species_fit$arrow_factor * 1,y=0,yend=PC2* species_fit$arrow_factor * 1),
                 arrow = arrow(length = unit(0.2, "cm")),colour="grey",alpha = 0.5) +
    geom_text_repel(data=species_fit$df,aes(x=PC1* species_fit$arrow_factor * 1,y=PC2* species_fit$arrow_factor * 1,label=name),
                    size=3,color = "grey",min.segment.length = 0) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) + 
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  ####NMDS -----------------------------------------------------------------------
  #####bray_curtis
  sol <- metaMDS(bray_curtis)
  stress  <- as.numeric(sprintf("%.2f",sol$stress))
  title = sprintf("Non-metric Multidimensional Scaling (%s) (%s)",name,"bray curtis distance")
  pc1 <- "NMDS 1"
  pc2 <- "NMDS 2"
  nmds <- data.frame(PC1=sol$points[,1],PC2=sol$points[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_bray_curtis_R2)~"P = "~.(adonis_bray_curtis_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_bray_curtis_P))
  line_3 = bquote("Stress:"~.(stress))
  nonmetric_label = list(line_1,line_2,line_3)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(nmds$PC2),0.90*max(nmds$PC2),0.80*max(nmds$PC2))
  
  p5 <- ggplot(nmds, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) +
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #####Artchison
  sol <- metaMDS(Artchison_dist)
  stress  <- as.numeric(sprintf("%.2f",sol$stress))
  title = sprintf("Non-metric Multidimensional Scaling (%s) (%s)",name,"Artchison distance")
  pc1 <- "NMDS 1"
  pc2 <- "NMDS 2"
  nmds <- data.frame(PC1=sol$points[,1],PC2=sol$points[,2],ind=alpha_diversity$sample,
                     Site=alpha_diversity$Site,Season = alpha_diversity$Season,Gradient = alpha_diversity$Gradient)
  
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_Artchison_R2)~"P = "~.(adonis_Artchison_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_Artchison_P))
  line_3 = bquote("Stress:"~.(stress))
  nonmetric_label = list(line_1,line_2,line_3)
  coord_x <- Inf#max(pca$PC1)
  coord_y <- c(max(nmds$PC2),0.85*max(nmds$PC2),0.70*max(nmds$PC2))
  
  p6 <- ggplot(nmds, aes(x=PC1,y=PC2))+
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    labs(x=pc1,y=pc2) +
    scale_shape_manual(values = c(17,16))+
    scale_colour_manual(values = cbPalette) +
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #final_plot <- (p1 + p2)/(p3 + p4) + plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A')
  return_list <- list(P_Shannon = p1,
                      P_Simpson = p2,
                      #P_pca = pca_plot,
                      P_PCoA_bray = p3,
                      P_PCoA_Artchison = p4,
                      P_NMDS_bray = p5,
                      P_NMDS_Artchison = p6
  )
  return(return_list)
} 
### Main part -----------------------------------------------------------------

COG_list_assembly_prok <- Diversity_2(COG_df_raw,COG,"COG_assembly_prok",taxonomy_list)
COG_path_list_assembly_prok <- Diversity_1(COG_pathway_df_raw,COG_path,"COG_path_assembly_prok",taxonomy_list)

pdf("~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiome_v3/Figures_v3/COG.pdf",
    width = 18,height = 6)
COG_list_assembly_prok$P_Shannon + 
  COG_list_assembly_prok$P_Simpson +
  COG_list_assembly_prok$P_PCoA_bray +
  plot_layout(guides = 'collect',nrow = 1)  &
  theme_prism(base_size = 16) &
  theme(legend.position='right')
dev.off()

return_list <- list(
                    COG_list_assembly_prok = COG_list_assembly_prok,
                    COG_path_list_assembly_prok = COG_path_list_assembly_prok
                    )
saveRDS(return_list, file="~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/Functional_alpha_beta_diversity_envfit_manuscript_list_2.rds")

return(return_list)
}
Functional_alpha_beta_diversity_list <- Functional_alpha_beta_diversity()
saveRDS(Functional_alpha_beta_diversity_list, file="~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Functional_alpha_beta_diversity_envfit_list.rds")








