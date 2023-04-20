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
  data_path <- "~/Desktop/marine_sediment/data/"
  cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02","CI" ="#e7298a","SW" = "#b10c0c","CDA"= "#ab37c9","BI" = "#e96e0f","TPC" = "#019a5b")
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  ### Setting color parameters ---------------------------------------------------
  Alpha_diversity <- function(data,name){
    # data <- paste0(data_path,"phylum.otu")
    # name <- "Phylum"
    taxon <- read_tsv(data) %>%
      dplyr::select(.,-PC_S_51,-PC_W_33) %>%
      column_to_rownames(var = "SampleID") %>%
      as.matrix() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      group_by(rowname) %>%
      rowwise() %>%
      mutate(Richness = sum(c_across(everything()) != 0)) %>%
      dplyr::select(rowname,Richness)
    alpha_diversity <- taxon %>%
      dplyr::rename(SampleID = rowname) %>%
      mutate(tmp = SampleID) %>%
      separate(tmp, c("Site", "Season","others")) %>%
      mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
      dplyr::select(.,-others) %>%
      as_tibble()
    gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                       Gradient = c("High","High","High","Medium","Medium","Low","Low"))
    alpha_diversity <- alpha_diversity %>%
      left_join(gradient)
    # The palette with grey:
    cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
    #### Figure 1
    kruskal.p <- kruskal.test(Richness~Gradient,data = alpha_diversity) %>% tidy() %>%
      .$p.value %>% round(2)
    kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
    tool2 = alpha_diversity %>%
      wilcox_test(Richness ~ Gradient) %>% 
      adjust_pvalue(method = 'fdr') %>%
      add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                       symbols = c("****", "***", "**", "*", "ns")
      )
    alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
    alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
    p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Richness",outlier.shape=NA) +
      ylab(paste("Richness","(",name,")",sep = "")) +
      geom_point(data = alpha_diversity, aes(x = Gradient, y = Richness, shape=Season,color = Site), position = position_jitter(width = 0.2))+
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = cbPalette) + xlab("")+
      annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
      stat_pvalue_manual(tool2, label = "p.adj.signif", 
                         hide.ns = TRUE,
                         size = 5,
                         y.position = seq(max(alpha_diversity$Richness),
                                          max(alpha_diversity$Richness)+0.05*length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]),
                                          length.out=length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]))
      ) +
      theme_prism(base_size = 10)
    return_result <- list(taxon = taxon,Shannon=p1)
    return(return_result)
  }
  Alpha_diversity_1 <- function(data,name){
    # data <- COG_path_assembly_prok
    # name <- "COG"
    taxon <- data %>%
      as.matrix() %>%
      as.data.frame() %>%
      mutate_all(funs(replace_na(.,0))) %>%
      rownames_to_column() %>%
      group_by(rowname) %>%
      rowwise() %>%
      mutate(Richness = sum(c_across(everything()) != 0)) %>%
      dplyr::select(rowname,Richness)
    alpha_diversity <- taxon %>%
      dplyr::rename(SampleID = rowname) %>%
      mutate(tmp = SampleID) %>%
      separate(tmp, c("Site", "Season","others")) %>%
      mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
      dplyr::select(.,-others) %>%
      as_tibble()
    gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                       Gradient = c("High","High","High","Medium","Medium","Low","Low"))
    alpha_diversity <- alpha_diversity %>%
      left_join(gradient)
    # The palette with grey:
    cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
    #### Figure 1
    kruskal.p <- kruskal.test(Richness~Gradient,data = alpha_diversity) %>% tidy() %>%
      .$p.value %>% round(2)
    kruskal.label <- paste("Kruskal-Wallis H test P value:",kruskal.p)
    tool2 = alpha_diversity %>%
      wilcox_test(Richness ~ Gradient) %>% 
      adjust_pvalue(method = 'fdr') %>%
      add_significance("p.adj",cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                       symbols = c("****", "***", "**", "*", "ns")
      )
    alpha_diversity$Site <- factor(alpha_diversity$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
    alpha_diversity$Gradient <- factor(alpha_diversity$Gradient,levels = c("Low","Medium","High"))
    p1 <- ggboxplot(alpha_diversity, x = "Gradient", y = "Richness",outlier.shape=NA) +
      ylab(paste("Richness","(",name,")",sep = "")) +
      geom_point(data = alpha_diversity, aes(x = Gradient, y = Richness, shape=Season,color = Site), position = position_jitter(width = 0.2))+
      scale_shape_manual(values = c(17,16)) +
      scale_color_manual(values = cbPalette) + xlab("")+
      annotate("text", x = Inf, y = Inf, label= kruskal.label,hjust = 1, vjust = 1,size = 4) + 
      stat_pvalue_manual(tool2, label = "p.adj.signif", 
                         hide.ns = TRUE,
                         size = 5,
                         y.position = seq(max(alpha_diversity$Richness),
                                          max(alpha_diversity$Richness)+0.05*length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]),
                                          length.out=length(tool2$p.adj.signif[tool2$p.adj.signif != "ns"]))
      ) +
      theme_prism(base_size = 10)
    return_result <- list(taxon = taxon,Shannon=p1)
    return(return_result)
  }
  #### Main part -----------------------------------------------------------------
  #superkingdom_list <- Alpha_diversity(paste0(data_path,"superkingdom.otu"),"Superkingdom")
  phylum_list <- Alpha_diversity(paste0(data_path,"phylum.otu"),"Phylum")
  class_list <- Alpha_diversity(paste0(data_path,"class.otu"),"Class")
  order_list <- Alpha_diversity(paste0(data_path,"order.otu"),"Order")
  family_list <- Alpha_diversity(paste0(data_path,"family.otu"),"Family")
  genus_list <- Alpha_diversity(paste0(data_path,"genus.otu"),"Genus")
  species_list <- Alpha_diversity(paste0(data_path,"species.otu"),"Species")
  
  ### Load raw data --------------------------------------------------------------
  # Assembly prok
  COG_file_assembly_prok <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_cog.tsv"
  COG_path_file_assembly_prok <- "~/Desktop/marine_sediment/RDS_files/functional/functional_assembly_metaspades/Assembly_prok_flt_v1/COG_pathway_stat.flt.txt"
  ### Read File ------------------------------------------------------------------
  COG_assembly_prok <- read_tsv(COG_file_assembly_prok) %>%
    column_to_rownames("Category") %>%
    t()
  COG_path_from_web <- read_tsv("~/Desktop/marine_sediment/results_new/functional_analyses/COG_pathway_from_web.tsv")
  COG_path_assembly_prok <- read_tsv(COG_path_file_assembly_prok) %>%
    filter(COG_pathway %in% COG_path_from_web$Pathway) %>%
    column_to_rownames("COG_pathway") %>%
    t()
  COG_list <- Alpha_diversity_1(COG_assembly_prok,"COG")
  #COG_path_list <- Alpha_diversity_1(COG_path_assembly_prok,"COG pathway")
  ### Small animal  ------------------------------------------------------------
  Animal_OTU_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/OTU_2019.csv")
  colnames(Animal_OTU_file)[1] <- "ID"
  Animal_OTU <- Animal_OTU_file %>%
    column_to_rownames("ID") %>%
    type_convert()
  Animal_OTU[Animal_OTU > 0] <- 1
  Animal_OTU <- Animal_OTU %>%
    rownames_to_column("ID")
  Animal_SAMP <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/SAMP_Linlin.csv")
  colnames(Animal_SAMP)[1] <- "ID"
  dim(Animal_SAMP)
  Animal_SAMP <- Animal_SAMP %>%
    dplyr::select(.,-ARMS_ID,-Ret_ID)
  Animal_TAXA_file <- read_csv("~/Desktop/marine_sediment/data/Shelby_data/Taxa_2019.csv")
  Animal_TAXA <- Animal_TAXA_file %>%
    filter(SourceID %in% "95_Match") %>%
    dplyr::select(.,-`...1`,-indtype,-SourceID) %>%
    rename(ID = Seq) %>%
    relocate(ID,.before = Kingdom)
  #colnames(Animal_TAXA)[1] <- "ID"
  dim(Animal_TAXA)
  Animal_TAXA_df <- Animal_TAXA %>%
    unite(label,Kingdom:Species)
  Animal_OTU_df <- Animal_OTU %>%
    left_join(Animal_SAMP) %>%
    relocate(c(Station_ID,Month_Ret), .after = ID) %>%
    mutate(Month_Ret = str_remove_all(Month_Ret,"[INUM]")) %>%
    mutate(ID = str_remove_all(ID,"ARMS_")) %>%
    unite("SampleID",c("Station_ID","Month_Ret","ID"))
  Animal_otu <- Animal_OTU_df %>%
    #filter(SampleID %in% rownames(COG_path)) %>%
    column_to_rownames("SampleID") 
  Animal_species_richness <- Animal_OTU_df %>%
    #filter(SampleID %in% rownames(COG_path)) %>%
    column_to_rownames("SampleID") %>% 
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname %in% Animal_TAXA$ID) %>%
    #left_join(Animal_TAXA,by = c("rowname"= "ID")) %>%
    #filter(!is.na(Species)) %>%
    #group_by(Species) %>%
    #summarise_if(is.double,max) %>%
    column_to_rownames() %>%
    t()  %>%
    as.data.frame() %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    rownames_to_column("SampleID") %>%
    mutate(tmp = SampleID) %>%
    dplyr::select(SampleID,tmp,total) %>%
    separate(tmp,c("Site","Season","tmp")) %>%
    dplyr::select(.,-tmp) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    left_join(gradient) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC"))) %>%
    mutate(Gradient = factor(Gradient, levels = c("High","Medium","Low")))
  Animal_species_richness <- Animal_species_richness %>%
    filter(SampleID %in% species_list$taxon$rowname) %>%
    arrange(SampleID)
  ## Plot ----------------------------------------------------------------------
  Spe_richness <- species_list$taxon %>%
    arrange(rowname) %>%
    column_to_rownames()
  Genus_richness <- genus_list$taxon %>%
    arrange(rowname) %>%
    column_to_rownames()
  Cog_richness <- COG_list$taxon %>%
    arrange(rowname) %>%
    column_to_rownames()
  Plot_df <- data.frame(Species = Spe_richness$Richness,
                        Genus = Genus_richness$Richness,
                        Cog = Cog_richness$Richness,
                        Animal = Animal_species_richness$total,
                        tmp = rownames(Spe_richness)) %>%
    `rownames<-`(rownames(Spe_richness)) %>%
    separate(tmp,c("Site","Season","tmp")) %>%
    dplyr::select(.,-tmp) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    mutate(Site = factor(Site,levels = c("SSW","PC","CI","SW","CDA","BI","TPC")))
  P1 <- ggplot(Plot_df,aes(x = Species, y = Cog)) +
    geom_point(aes(color = Site,shape = Season)) +
    geom_rug() +
    geom_smooth(method = lm) +
    ggpubr::stat_cor(size = 12/.pt) +
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) +
    ylab("COG Richness") +
    xlab("Species Richness") +
    theme_prism(base_size = 16)
  P2 <- ggplot(Plot_df,aes(x = Species, y = Animal)) +
    geom_point(aes(color = Site,shape = Season)) +
    geom_rug() +
    geom_smooth(method = lm) +
    ggpubr::stat_cor(size = 12/.pt) +
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) +
    ylab("Benthos Richness") +
    xlab("Species Richness") +
    theme_prism(base_size = 16)
  P3 <- ggplot(Plot_df,aes(x = Genus, y = Animal)) +
    geom_point(aes(color = Site,shape = Season)) +
    geom_rug() +
    geom_smooth(method = lm) +
    ggpubr::stat_cor(size = 12/.pt) +
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) +
    ylab("Benthos Richness") +
    xlab("Genus Richness") +
    theme_prism(base_size = 16)
  P4 <- ggplot(Plot_df,aes(x = Cog, y = Animal)) +
    geom_point(aes(color = Site,shape = Season)) +
    geom_rug() +
    geom_smooth(method = lm) +
    ggpubr::stat_cor(size = 12/.pt) +
    scale_shape_manual(values = c(17,16)) +
    scale_color_manual(values = cbPalette) +
    ylab("Benthos Richness") +
    xlab("COG Richness") +
    theme_prism(base_size = 16)
  p <- P2 + P3 + P4  + 
    plot_layout(ncol = 3) +
    plot_layout(guides = 'collect')  &
    theme(legend.position='bottom') 
  
  cor.test(Spe_richness$Richness,Cog_richness$Richness)
  plot(Spe_richness$Richness,Cog_richness$Richness,
       xlab = "Species Richness", ylab = "COG Richness",
       col = rgb(0.5,0.5,0.8),cex = 1,pch = 19)
  lmodel <- lm(Cog_richness$Richness~Spe_richness$Richness)
  abline(lmodel,lwd = 2)
  pdf("~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/Compare_richness_otu_95match.pdf",width=12, height= 4)
  p
  dev.off()
