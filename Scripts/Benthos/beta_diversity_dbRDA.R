# Title: Bata diversity benthos
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

### Loading required package ---------------------------------------------------
library(vegan)
library(ggplot2)
library(ggnewscale)
library(reshape)
library(dplyr)
library(vegan)
library(robCompositions)
library(devtools)
library(amplicon)
library(ape)
library(patchwork)
library(tibble)
library(grid)
library(RColorBrewer)
library(packfor)
library(ggrepel)
library(tidyverse)
library(ggprism)
library(phyloseq)
### Function Beta_diversity ----------------------------------------------------
Beta_diversity <- function(data)
{
  ## Tidy data -----------------------------------------------------------------
  # data <- Animal_otu
  
  taxon <- data
  
  taxon1 <- taxon %>%
    column_to_rownames(var = "Feature") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
  df <- taxon1 %>%
    as.data.frame() %>%
    rownames_to_column(var = "SampleID") %>%
    mutate(tmp = SampleID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble()
  ## beta-diversity ------------------------------------------------------------
  ## Bray_Curtis
  taxon2 <- as.matrix(taxon1)
  jaccard_dist <- vegdist(taxon2, method = "jaccard")
  jaccard_dist <- as.matrix(jaccard_dist)
  # PERMANOVA + betadisper + anosim analysis
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","High","Low","Low","Low"))
  df <- df %>%
    left_join(gradient)
  all(rownames(jaccard_dist) == df$SampleID)
  all(rownames(taxon2) == rownames(jaccard_dist))
  ###PERMANOVA
  adonis_jaccard_dist = adonis2(jaccard_dist~Gradient, df, permutations = 999)
  adonis_jaccard_R2 = round(adonis_jaccard_dist$R2[1],2)
  adonis_jaccard_P = round(adonis_jaccard_dist$`Pr(>F)`[1],4)
  
  adonis_jaccard_dist_1 = adonis2(jaccard_dist~Season, df, permutations = 999)
  adonis_jaccards_R2_1 = round(adonis_jaccard_dist_1$R2[1],2)
  adonis_jaccard_P_1 = round(adonis_jaccard_dist_1$`Pr(>F)`[1],4)
  ###betadisper
  beta_jaccard_dist <- betadisper(as.dist(jaccard_dist), df$Gradient) %>% permutest()
  beta_jaccard_P <- round(as.data.frame(beta_jaccard_dist$tab)[1,]$`Pr(>F)`,4)
  ## color palette -------------------------------------------------------------
  # The palette with grey:
  envPalette <- c("High" = "#D55E00","Low" = "#0072B2")
  cbPalette <- c("SSW" = "#0c75d1", "PC" = "#e6ab02", "CI" = "#e7298a", 
                    "SW" = "#b10c0c", "CDA" = "#ab37c9", "BI" = "#e96e0f",
                    "TPC" = "#019a5b"
  )
  ## Pollution data ------------------------------------------------------------
  file1 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ms5.csv")
  file2 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ns6.csv")
  file3 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ps6.csv")
  file4 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ss1-ss3-ss5.csv")
  file5 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5y-ts2.csv")
  file6 <- read_csv("~/Desktop/marine_sediment/2019 Water and Sediment Data/marine_pollution/marine_sedi_quality-5yrs-ms8.csv")
  file_sediment <- rbind(file1,file2,file3,file4,file5,file6)
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
    unite("group", StationID:season)
  df_form <- data.frame(sample = df$SampleID,
                        Site = df$Site,
                        Season = df$Season,
                        group = paste(df$Site,df$Season,sep = "_"))
  merge_tab <- df_form %>%
    left_join(Score_sediment)
  rownames(merge_tab) <- merge_tab$sample
  merge_tab_1 = merge_tab[,!(names(merge_tab) %in% c("group","sample","Site","Season"))]
  ## PcoA ----------------------------------------------------------------------
  ## Bray Curtis
  PCOA <- pcoa(jaccard_dist, correction="none", rn=NULL)
  PCOA1 <- data.frame(PC1 = PCOA$vectors[,1],PC2 = PCOA$vectors[,2])
  fit <- envfit(PCOA1,merge_tab_1,perm=999)
  fit.df<-as.data.frame(fit$vectors$arrows*sqrt(fit$vectors$r))
  fit.df$pollution<-rownames(fit.df)
  fit.df$Rsquare <- fit$vectors$r
  fit.df$pvals <- fit$vectors$pvals
  fit.df.1 <- fit.df %>% 
    mutate(BH = p.adjust(pvals)) %>%
    filter(pvals < 0.05)
  fit.df.2 <- fit.df %>% 
    filter(pvals >= 0.05)
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=df$SampleID,
                     Site=df$Site,Season = df$Season)
  pcoa <- pcoa %>%
    left_join(gradient)
  pcoa$Site <- factor(pcoa$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  pcoa$Gradient <- factor(pcoa$Gradient,levels = c("High","Low"))
  line_1 = bquote("PERMANOVA :" ~ R^2 ~"= "~.(adonis_jaccard_R2)~"P = "~.(adonis_jaccard_P))
  line_2 = bquote("PERMANOVA season:" ~ R^2 ~"= "~.(adonis_jaccards_R2_1)~"P = "~.(adonis_jaccard_P_1))
  #line_3 = bquote("PERMDISP2:"~"P = "~.(beta_bray_curtis_P))
  #nonmetric_label = list(line_1,line_2,line_3)
  nonmetric_label = list(line_1,line_2)
  coord_x <- max(pcoa$PC1)
  #coord_y <- c(max(pcoa$PC2),0.9*max(pcoa$PC2),0.8*max(pcoa$PC2))
  coord_y <- c(max(pcoa$PC2),0.9*max(pcoa$PC2))
  arrow_factor <- ordiArrowMul(fit) *0.2
  #envPalette <- c("#FB8072","#FDB462","#8DD3C7")
  p1 <- ggplot(data = pcoa, aes(x=PC1,y=PC2))+#,color=pcoa$Site,shape = pcoa$Season))+
    geom_point(data = pcoa, aes(color=Site,shape = Season),size=2,position = "jitter") + 
    labs(x=pc1,y=pc2) + 
    geom_segment(data=fit.df.1,aes(x=0,xend=PC1* arrow_factor,y=0,yend=PC2* arrow_factor,color = col),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    ggrepel::geom_text_repel(data=fit.df.1,aes(x=PC1* arrow_factor * 1.1,y=PC2* arrow_factor* 1.1,label=pollution),
                             size=3,colour="#EEAD0E")+
    # geom_segment(data=fit.df.2,aes(x=0,xend=PC1* arrow_factor,y=0,yend=PC2* arrow_factor,color = col),
    #              arrow = arrow(length = unit(0.2, "cm")),colour="#6495ED",alpha = 0.5) +
    # geom_text(data=fit.df.2,aes(x=PC1* arrow_factor * 1.1,y=PC2* arrow_factor* 1.1,label=pollution),size=3,colour="#6495ED")+
    geom_hline(yintercept = 0,linetype=2,color="grey")+geom_vline(xintercept = 0,linetype=2,color="grey")+
    #geom_path(data=conf.rgn)+
    scale_colour_manual(values = cbPalette) + 
    scale_shape_manual(values = c(17,16))+
    labs(title = title)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1) +
    theme_bw()+
    #stat_ellipse(data = pcoa,aes(group = Gradient,fill=Gradient), type = "norm",level = 0.68)+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  #p1
  
  ##dbRDA-----
  merge_tab_2 <- merge_tab_1
  R_square_table <- data.frame()
  rm(formula_cum)
  for (i in colnames(merge_tab_1)) {
    i1 <- i
    i <- paste0("`",i,"`")
    formula <- as.formula(paste0("taxon2~",i))
    set.seed(123)
    uni_dbrda <- capscale(formula, merge_tab_2, distance = 'jaccard')
    uni_i <- max(RsquareAdj(uni_dbrda)$r.squared,0)
    uni_p <- anova(uni_dbrda)$`Pr(>F)`[1]
    if (!exists("formula_cum")) {
      formula_cum <- formula
    }else{
      formula_cum <- as.formula(paste0(Reduce(paste, deparse(formula_cum)),"+",i))
    }
    set.seed(123)
    db_rda <- capscale(formula_cum,
                       data = merge_tab_2,dist = "jaccard")
    cum_i <- RsquareAdj(db_rda)$r.squared
    re_tb <- data.frame(name = i1, uni = uni_i,unip  = uni_p, cum = cum_i)
    R_square_table <- rbind(R_square_table,re_tb)
  }
  pollut_category <- read_tsv("~/Desktop/marine_sediment/renew_pollutionbin/sediment_pollute.tsv")
  
  R_square_table_cate <- R_square_table %>%
    left_join(pollut_category) %>%
    mutate(BHP = p.adjust(unip,method = "BH"))
  R_sign_table <- R_square_table_cate %>% 
    filter(BHP < 0.05)

  set.seed(123)
  db_rda <- capscale(taxon2 ~ .,
                     data = merge_tab_2,dist = "jaccard")
  db_rda_r <- RsquareAdj(db_rda)$r.squared
  db_rda_r.adjusted <- RsquareAdj(db_rda)$adj.r.squared
  anova(db_rda)
  db_rda_p <- anova(db_rda)$`Pr(>F)`
  
  #forward stepwise selection, ordiR2step() permutations: 999
  set.seed(123)
  # db_rda_backward_pr <- ordiR2step(capscale(taxon~1, merge_tab_2, distance = 'bray'), scope = formula(db_rda),
  #                              R2scope = F, direction = 'backward', permutations = 999)
  # db_rda_both_pr <- ordiR2step(capscale(taxon~1, merge_tab_2, distance = 'bray'), scope = formula(db_rda),
  #                                 R2scope = F, direction = 'both', permutations = 999)
  db_rda_forward_pr <- ordiR2step(capscale(taxon2~1, merge_tab_2, distance = 'jaccard'), scope = formula(db_rda),
                                  R2scope = T, direction = 'forward', permutations = 999)
  
  #db_rda_both_pr <- ordiR2step(capscale(taxon~1, merge_tab_2, distance = 'jaccard', add = TRUE), scope = formula(db_rda), direction = 'both', permutations = 999)
  db_rda_forward_pr_R <- RsquareAdj(db_rda_forward_pr)$r.squared

  levels <- db_rda_forward_pr$anova %>% 
    as.data.frame() %>% 
    rownames() %>%
    str_remove("\\+ ") %>%
    str_remove_all('`') %>%
    head(-1)
  
  re_list <- list(
    R_square_table_cate = R_square_table_cate,
    forward_select = levels,
    db_rda_forward_pr_R = db_rda_forward_pr_R,
    db_rda_r = db_rda_r,
    db_rda_forward_pr = db_rda_forward_pr
  )
  #saveRDS(re_list,"~/Desktop/marine_sediment/results/MEGAN6/dbRDA.rds")
  #Permutation Tests of RDA Results
  db_rda_forward_pr_test <- anova(db_rda_forward_pr, permutations = 999)
  
  #More details
  #summary(db_rda_forward_pr, scaling = 1)
  
  # Loadings of environment parameter and microbiomes, scaling = 1
  db_rda_forward_pr.scaling1 <- summary(db_rda_forward_pr, scaling = 1)
  db_rda_forward_pr.site <- data.frame(db_rda_forward_pr.scaling1$sites)[1:2]
  db_rda_forward_pr.env <- data.frame(db_rda_forward_pr.scaling1$biplot)[1:2]
  #db_rda_forward_pr.spe <- data.frame(db_rda_forward_pr.scaling1$species)[1:2] * 0.2
  #db_rda_forward_pr.spe <- db_rda_forward_pr.spe[colnames(taxon_top_10),]
  #add group name
  db_rda_forward_pr.env$name <- rownames(db_rda_forward_pr.env) 
  db_rda_forward_pr.env <- db_rda_forward_pr.env %>%
    mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
    mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
    mutate(name = str_remove_all(name,"`"))
  #db_rda_forward_pr.spe$name <- rownames(db_rda_forward_pr.spe)
  db_rda_forward_pr.site$name <- rownames(db_rda_forward_pr.site)
  db_rda_forward_pr.site <- db_rda_forward_pr.site %>%
    mutate(tmp = name) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    left_join(gradient) %>%
    dplyr::select(.,-others)
  
  #Adjust R^2
  #exp_adj <- RsquareAdj(db_rda_forward_pr)$r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
  exp_adj <- db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
  rda1_exp <- paste('dbRDA1:', round(exp_adj[1]*100, 2), '%')
  rda2_exp <- paste('dbRDA2:', round(exp_adj[2]*100, 2), '%')
  #ggplot2
  library(ggplot2)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_jaccard_R2)~"P = "~.(adonis_jaccard_P))
  line_2 = bquote("PERMANOVA season:" ~ R^2 ~"= "~.(adonis_jaccards_R2_1)~"P = "~.(adonis_jaccard_P_1))
  #line_3 = bquote("PERMDISP2:"~"P = "~.(beta_bray_curtis_P))
  #line_4 = bquote("dbRDA anova:"~"F"~"= "~.(db_rda_forward_pr_test$F[1])~"P = "~.(db_rda_forward_pr_test$`Pr(>F)`[1]))
  line_4 = bquote("Total explained variation:"~ R^2 ~"= "~.(paste(round(db_rda_r*100, 2), '%')))
  line_5 = bquote("Constrained explained variation:"~ R^2 ~"= "~.(paste(round(db_rda_forward_pr_R*100, 2), '%')))
  nonmetric_label = list(line_1,line_2,line_4,line_5)
  #nonmetric_label = list(line_1,line_5)
  #coord_x <- 1.0 *max(db_rda_forward_pr.site$CAP1)
  coord_x <- Inf
  # coord_y <- c(1.2 * max(db_rda_forward_pr.site$CAP2),
  #              1.05 * max(db_rda_forward_pr.site$CAP2),
  #              0.9 * max(db_rda_forward_pr.site$CAP2),
  #              0.8 * max(db_rda_forward_pr.site$CAP2),
  #              0.7 * max(db_rda_forward_pr.site$CAP2))
  coord_y <- c(1.5 * max(db_rda_forward_pr.site$CAP2),
               1.38 * max(db_rda_forward_pr.site$CAP2),
               1.26 * max(db_rda_forward_pr.site$CAP2),
               1.14 * max(db_rda_forward_pr.site$CAP2)
               )
  title = sprintf("Stepwise dbRDA")
  p2 <- ggplot(db_rda_forward_pr.site, aes(CAP1, CAP2)) +
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    scale_color_manual(name = "Site",values = cbPalette) +
    new_scale_color() +
    #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
    labs(x = rda1_exp, y = rda2_exp) +
    geom_hline(yintercept = 0,linetype=2,color="grey")+geom_vline(xintercept = 0,linetype=2,color="grey")+
    #geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_shape_manual(values = c(17,16))+
    labs(title = title)+
    geom_segment(data = db_rda_forward_pr.env, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
                 arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'black') +
    # geom_segment(data = db_rda_forward_pr.spe, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
    #              arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = '#6495ED') +
    #geom_text(data = db_rda_forward_pr.spe, aes(CAP1 * 1.1, CAP2 * 1.1, label = name), color = 'grey', size = 3) +
    # geom_text_repel(data=db_rda_forward_pr.spe,aes(x = CAP1*1, 
    #                                                y = CAP2 * 1, 
    #                                                label = name),
    #                 size=5,min.segment.length = unit(0, 'lines'),
    #                 color = "#6495ED") +
    geom_text_repel(data=db_rda_forward_pr.env,aes(x = CAP1*1, 
                                                   y = CAP2 * 1, 
                                                   label = name),
                    size=5,min.segment.length = unit(0, 'lines'),
                    color = "black") +
    #theme_bw() +
    theme_prism(base_size = 16,border = T)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),
             hjust = 1, vjust = 1,size = 16/.pt) +
    #geom_text(data = db_rda_forward_pr.env, aes(CAP1 * 1.1, CAP2 * 1.1, label = name), color = 'blue', size = 3) +
    stat_ellipse(aes(Site = Gradient,color=Gradient),
                 level=0.95,alpha=1,
                 type = "norm")+
    scale_color_manual(name = "Site",values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #p2
  # Unconstrained variation showed in PCoA ------
  res_forward_pr.site <- data.frame(db_rda_forward_pr$CA$u)[1:2]
  res_forward_pr.site$name <- rownames(res_forward_pr.site)
  res_forward_pr.site <- res_forward_pr.site %>%
    mutate(tmp = name) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    left_join(gradient) %>%
    dplyr::select(.,-others)

  rda1_exp <- paste('PCOA1:', round((db_rda$CA$eig/sum(db_rda$CA$eig))[1]*100, 2), '%')
  rda2_exp <- paste('PCOA2:', round((db_rda$CA$eig/sum(db_rda$CA$eig))[2]*100, 2), '%')
  title = sprintf("Unconstrained ordination ")
  p3 <- ggplot(db_rda_forward_pr.site, aes(CAP1, CAP2)) +
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    scale_color_manual(values = cbPalette) +
    labs(x = rda1_exp, y = rda2_exp) +
    geom_hline(yintercept = 0,linetype=2,color="grey")+geom_vline(xintercept = 0,linetype=2,color="grey")+
    scale_shape_manual(values = c(17,16))+
    theme_bw() +
    labs(title = title)+
    stat_ellipse(aes(group = Gradient,fill=Gradient),
                 geom="polygon",level=0.68,alpha=0.1,
                 type = "norm")+
    scale_fill_manual(values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  #p3
  # Calculate environment parameter contribution. ---------
  levels <- db_rda_forward_pr$anova %>% 
    as.data.frame() %>% 
    rownames() %>%
    str_remove("\\+ ") %>%
    str_remove_all('`') %>%
    head(-1)
  Contribut_df <- merge_tab_1 %>%
    relocate(levels, .before = 1)
  R_square_table <- data.frame()
  rm(formula_cum)
  for (i in colnames(Contribut_df)) {
    i1 <- i
    i <- paste0("`",i,"`")
    formula <- as.formula(paste0("taxon2~",i))
    set.seed(123)
    uni_dbrda <- capscale(formula, Contribut_df, distance = 'jaccard')
    uni_i <- max(RsquareAdj(uni_dbrda)$r.squared,0)
    uni_p <- anova(uni_dbrda)$`Pr(>F)`[1]
    if (!exists("formula_cum")) {
      formula_cum <- formula
    }else{
      formula_cum <- as.formula(paste0(Reduce(paste, deparse(formula_cum)),"+",i))
    }
    set.seed(123)
    db_rda <- capscale(formula_cum,
                       data = Contribut_df,dist = "jaccard")
    cum_i <- RsquareAdj(db_rda)$r.squared
    re_tb <- data.frame(name = i1, uni = uni_i,unip  = uni_p, cum = cum_i)
    R_square_table <- rbind(R_square_table,re_tb)
  }
  R_sign_table <- R_square_table %>% 
    mutate(BHP = p.adjust(unip,method = "BH")) %>% 
    filter(BHP < 0.05)
  return_result <- list(pcoa=p1,dbRDA=p2,pcoa_dbRDA=p3,
                        cumRdf = R_square_table,levels = levels,
                        R_sign_table = R_sign_table
                        )
  return(return_result)
}
### Function Barplot -----------------------------------------------------------
Barplot <- function(data,sign_v)
{
  # data <- Species_list$cumRd
  # sign_v <- Species_list$levels
  data <- data %>%
    mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
    mutate(name = str_remove_all(name," \\(% Total Solid\\)")) %>%
    mutate(name = factor(name,levels = rev(name))) %>%
    pivot_longer(2:4,names_to = "group",values_to = "value")
  
  sign_ID <- data.frame(name = sign_v) %>%
    mutate(name = str_remove_all(name," \\(μg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/kg\\)")) %>%
    mutate(name = str_remove_all(name," \\(mg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(μg/L\\)")) %>%
    mutate(name = str_remove_all(name," \\(%w/w\\)")) %>%
    mutate(name = str_remove_all(name," \\(% Total Solid\\)"))
  
  sign_R <- data %>%
    filter(name %in% tail(sign_ID,1),group %in% "cum") %>%
    .$value
  data <- data %>%
    filter(!group == "unip") 
  p3 <- ggplot(data) +
    geom_bar(stat = "identity",position = "dodge",aes(x=name,y=value,fill = group,group = group)) +
    scale_fill_manual(values = c("grey","black"),labels = c(bquote("Cumulative" ~ R^2),bquote("Univariate" ~ R^2))) +
    coord_flip() + geom_hline(yintercept = sign_R,linetype=2,color="darkred")+
    ylab("Effect size") + xlab("Covariates of microbiome composition") +
    #theme_bw() + 
    theme_prism(base_size = 16,border = T) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(legend.title = element_blank()) +
    theme(legend.position="bottom")
  return(p3)
}
### Main part -----------------------------------------------------------------
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
Animal_TAXA_file <- read_tsv("~/Desktop/marine_sediment/data/Shelby_data/TAXA_2019_onlysediment.tsv")
Animal_TAXA <- Animal_TAXA_file %>%
  dplyr::select(.,-ID,-indtype,-SourceID) %>%
  dplyr::rename(ID = Seq) %>%
  relocate(ID,.before = Kingdom) %>%
  mutate(OTU = Species) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Genus),paste0(Genus," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Family),paste0(Family," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Order),paste0(Order," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Class),paste0(Class," OTU",ID),OTU)) %>%
  mutate(OTU = if_else(is.na(OTU) & !is.na(Phylum),paste0(Phylum," OTU",ID),OTU))
#colnames(Animal_TAXA)[1] <- "ID"
dim(Animal_TAXA)
Animal_OTU_df <- Animal_OTU %>%
  left_join(Animal_SAMP) %>%
  relocate(c(Station_ID,Month_Ret), .after = ID) %>%
  mutate(Month_Ret = str_remove_all(Month_Ret,"[INUM]")) %>%
  mutate(ID = str_remove_all(ID,"ARMS_")) %>%
  unite("SampleID",c("Station_ID","Month_Ret","ID"))

Animal_otu <- Animal_OTU_df %>%
  #filter(SampleID %in% rownames(COG_path)) %>%
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  filter(Feature %in% Animal_TAXA$ID) %>%
  as.data.frame()

rm(formula_cum)
OTU_list <- Beta_diversity(Animal_otu)
#### Species level -------------------------------------------------------------
bar_plot <- Barplot(OTU_list$cumRdf,OTU_list$levels)
OTU_list$dbRDA_withbarplot <- OTU_list$dbRDA + bar_plot +
  plot_layout(design = "AAAB") +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 

pdf("~/Desktop/marine_sediment/results/Benthos/dbRDA_OTU1.pdf",width = 18,height = 10)
#beta_diversity_dbRDA_list$Species$dbRDA_withbarplot
OTU_list$dbRDA_withbarplot
dev.off()

saveRDS(OTU_list,"~/Desktop/marine_sediment/results/Benthos/dbRDA_2019_OTU.rds")

