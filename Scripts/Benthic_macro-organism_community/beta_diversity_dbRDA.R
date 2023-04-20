beta_diversity_dbRDA <- function()
  {
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
# cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
# "#D55E00", "#CC79A7")
### Function Beta_diversity ----------------------------------------------------
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
    filter(pvals < 0.05,pvals.BH < 0.1)
  fit.df.2 <- fit.df %>% 
    filter(pvals >= 0.05)
  arrow_factor <- ordiArrowMul(fit)
  return_list <- list(
    df = fit.df.1,
    tmp = fit.df,
    arrow_factor = arrow_factor
  )
}
Beta_diversity <- function(data,name)
{
  ## Tidy data -----------------------------------------------------------------
  data <- Animal_species
  name <- "Species"
  taxon <- data %>%
    column_to_rownames(var = "ID") #%>%
    # as.matrix() %>%
    # t() %>%
    # as.data.frame()
  df <- data %>%
    mutate(tmp = ID) %>%
    separate(tmp, c("Site", "Season","others")) %>%
    mutate(Season = str_replace_all(Season,c("S" = "Summer", "W" = "Winter"))) %>%
    dplyr::select(.,-others) %>%
    as_tibble()
  ## beta-diversity ------------------------------------------------------------
  ####bray_curtis
  taxon <- as.matrix(taxon)
  jaccard_dis <- vegdist(taxon, method = "jaccard")
  jaccard_dis <- as.matrix(jaccard_dis)
  # PERMANOVA + betadisper + anosim analysis
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  df <- df %>%
    left_join(gradient)
  ###PERMANOVA
  adonis_jaccard_dist = adonis2(jaccard_dis~Gradient, df, permutations = 999)
  adonis_jaccard_R2 = round(adonis_jaccard_dist$R2[1],2)
  adonis_jaccard_P = round(adonis_jaccard_dist$`Pr(>F)`[1],4)
  ###betadisper
  beta_jaccard_dist <- betadisper(as.dist(jaccard_dis), df$Gradient) %>% permutest()
  beta_jaccard_P <- round(as.data.frame(beta_jaccard_dist$tab)[1,]$`Pr(>F)`,4)
  ## color palette -------------------------------------------------------------
  # The palette with grey:
  cbPalette <- c("#0c75d1", "#e6ab02", "#e7298a", "#b10c0c", "#ab37c9", "#e96e0f", "#019a5b")
  envPalette <- c("#D55E00","#E69F00","#0072B2")
  ## Pollution data ------------------------------------------------------------
  pollution <- read_tsv("~/Desktop/marine_sediment/results/pollution/year2019_pollution.tsv")
  pollution_form <- pollution %>%
    relocate(c("Nitrate Nitrogen (mg/L)","Nitrite Nitrogen (mg/L)","Total Inorganic Nitrogen (mg/L)","Total Nitrogen (mg/L)"),
             .after = "season") %>%
    dplyr::select(.,-year) %>%
    unite("group", StationID:season)
  df_form <- data.frame(sample = rownames(df),
                        Site = df$Site,
                        Season = df$Season,
                        group = paste(df$Site,df$Season,sep = "_"))
  merge_tab <- merge(df_form,pollution_form)
  rownames(merge_tab) <- merge_tab$sample
  colnames(merge_tab) <- str_remove_all(colnames(merge_tab)," \\(μg/kg\\)") %>%
    str_remove_all(.," \\(mg/kg\\)") %>%
    str_remove_all(.," \\(mg/L\\)") %>%
    str_remove_all(.," \\(μg/L\\)") %>%
    str_remove_all(.," \\(%w/w\\)") %>%
    str_remove_all(.," \\(% Total Solid\\)")
  merge_tab_1 = merge_tab[,!(names(merge_tab) %in% c("group","sample","Site"))]
  merge_tab = merge_tab[,!(names(merge_tab) %in% c("group","sample","Site","Season"))]
  ## PcoA ----------------------------------------------------------------------
  #####jaccard
  PCOA <- pcoa(jaccard_dis, correction="none", rn=NULL)
  PCOA1 <- data.frame(PC1 = PCOA$vectors[,1],PC2 = PCOA$vectors[,2])
  fit_list <- func_envfit(PCOA1,merge_tab_1)
  fit.df.1 <- fit_list$df
  result <-PCOA$values[,"Relative_eig"]
  title = sprintf("Principal Coordinates Analysis (%s) (%s)",name,"jaccard distance")
  pc1 <- paste("PCoA1 (",round(result[1]*100,digits=2),"%)",sep="")
  pc2 <- paste("PCoA2 (",round(result[2]*100,digits=2),"%)",sep="")
  pcoa <- data.frame(PC1=PCOA$vectors[,1],PC2=PCOA$vectors[,2],ind=row.names(df),
                     Site=df$Site,Season = df$Season)
  gradient <- tibble(Site = c("SSW","PC","CI","SW","CDA","BI","TPC"),
                     Gradient = c("High","High","High","Medium","Medium","Low","Low"))
  pcoa <- pcoa %>%
    left_join(gradient)
  pcoa$Site <- factor(pcoa$Site,levels = c("SSW","PC", "CI", "SW", "CDA", "BI", "TPC"))
  pcoa$Gradient <- factor(pcoa$Gradient,levels = c("High","Medium","Low"))
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_jaccard_R2)~"P = "~.(adonis_jaccard_P))
  line_2 = bquote("PERMDISP2:"~"P = "~.(beta_jaccard_P))
  #line_3 = bquote("ANOSIM:"~"R"~"= "~.(anosim_bray_curtis_R)~"P = "~.(anosim_bray_curtis_P))
  nonmetric_label = list(line_1,line_2)
  coord_x <- max(pcoa$PC1)
  coord_y <- c(max(pcoa$PC2),0.85*max(pcoa$PC2))
  arrow_factor <- fit_list$arrow_factor *0.2
  #taxon_arrow_factor <- ordiArrowMul(taxon_fit) *0.25
  #envPalette <- c("#FB8072","#FDB462","#8DD3C7")
  p1 <- ggplot(data = pcoa, aes(x=PC1,y=PC2))+#,color=pcoa$Site,shape = pcoa$Season))+
    geom_point(data = pcoa, aes(color=Site,shape = Season),size=2,position = "jitter") + 
    labs(x=pc1,y=pc2) + 
    geom_segment(data=fit.df.1,aes(x=0,xend=PC1* arrow_factor,y=0,yend=PC2* arrow_factor),
                 arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
    geom_text_repel(data=fit.df.1,aes(x=PC1* arrow_factor ,y=PC2* arrow_factor,label=name),
                    size=3,colour="#EEAD0E",max.overlaps = 30,min.segment.length = 0)+
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
  if(dim(fit.df.1)[1] == 0){
    p1 <- ggplot(data = pcoa, aes(x=PC1,y=PC2))+#,color=pcoa$Site,shape = pcoa$Season))+
      geom_point(data = pcoa, aes(color=Site,shape = Season),size=2,position = "jitter") + 
      labs(x=pc1,y=pc2) + 
      # geom_segment(data=fit.df.1,aes(x=0,xend=PC1* arrow_factor,y=0,yend=PC2* arrow_factor),
      #              arrow = arrow(length = unit(0.2, "cm")),colour="#EEAD0E",alpha = 0.5) +
      # geom_text_repel(data=fit.df.1,aes(x=PC1* arrow_factor ,y=PC2* arrow_factor,label=name),
      #                 size=3,colour="#EEAD0E",max.overlaps = 30,min.segment.length = 0)+
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
  }
  p1
  
  ##dbRDA-----
  merge_tab_2 <- merge_tab_1 %>% `rownames<-`(rownames(taxon))
  #paste0("taxon ~ ",)
  set.seed(123)
  db_rda <- capscale(taxon ~ .,
                     data = merge_tab_2,dist = "jaccard")
  db_rda_r <- RsquareAdj(db_rda)$r.squared
  db_rda_r.adjusted <- RsquareAdj(db_rda)$adj.r.squared
  anova(db_rda)
  set.seed(123)
  #forward stepwise selection, ordiR2step() permutations: 999
  db_rda_forward_pr <- ordiR2step(capscale(taxon~1, merge_tab_2, distance = 'jaccard'), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
  #db_rda_both_pr <- ordiR2step(capscale(taxon~1, merge_tab_2, distance = 'bray', add = TRUE), scope = formula(db_rda), direction = 'both', permutations = 999)
  db_rda_forward_pr_Radj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared
  
  #Permutation Tests of RDA Results
  db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)
  
  #Permutation Tests of each axes
  db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)
  #db_rda_forward_pr_test_terms <- anova.cca(db_rda_forward_pr, by = 'terms', permutations = 999)
  #p value adjust
  db_rda_forward_pr_test_axis$`Pr(>F)` <- p.adjust(db_rda_forward_pr_test_axis$`Pr(>F)`, method = 'bonferroni')
  
  
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
  #db_rda_forward_pr.spe$name <- rownames(db_rda_forward_pr.spe)
  db_rda_forward_pr.site$name <- rownames(db_rda_forward_pr.site)
  db_rda_forward_pr.site$Site <- pcoa$Site
  db_rda_forward_pr.site$Season <- pcoa$Season
  db_rda_forward_pr.site$Gradient <- pcoa$Gradient
  
  #Adjust R^2
  exp_adj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
  rda1_exp <- paste('dbRDA1:', round(exp_adj[1]*100, 2), '%')
  rda2_exp <- paste('dbRDA2:', round(exp_adj[2]*100, 2), '%')
  #ggplot2
  library(ggplot2)
  line_1 = bquote("PERMANOVA:" ~ R^2 ~"= "~.(adonis_jaccard_R2)~"P = "~.(adonis_jaccard_P))
  #line_2 = bquote("PERMDISP2:"~"P = "~.(beta_jaccard_P))
  #line_3 = bquote("dbRDA anova:"~"F"~"= "~.(db_rda_forward_pr_test$F[1])~"P = "~.(db_rda_forward_pr_test$`Pr(>F)`[1]))
  line_4 = bquote("Constrained explained variation:"~ R^2 ~"= "~.(paste(round(db_rda_forward_pr_Radj*100, 2), '%')))
  #nonmetric_label = list(line_1,line_2,line_3,line_4)
  nonmetric_label = list(line_1,line_4)
  #coord_x <- 1.0 *max(db_rda_forward_pr.site$CAP1)
  coord_x <- Inf
  # coord_y <- c(1.0 * max(db_rda_forward_pr.site$CAP2),
  #              0.9 *max(db_rda_forward_pr.site$CAP2),
  #              0.8*max(db_rda_forward_pr.site$CAP2),
  #              0.7*max(db_rda_forward_pr.site$CAP2))
  coord_y <- c(1.0 * max(db_rda_forward_pr.site$CAP2),
               0.9 *max(db_rda_forward_pr.site$CAP2))
  title = sprintf("Stepwise dbRDA (%s) (%s)",name,"jaccard distance")
  p2 <- ggplot(db_rda_forward_pr.site, aes(CAP1, CAP2)) +
    geom_point(aes(color = Site,shape=Season),position = "jitter") +
    scale_color_manual(values = cbPalette) +
    new_scale_color() +
    #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
    labs(x = rda1_exp, y = rda2_exp) +
    geom_hline(yintercept = 0,linetype=2,color="grey")+geom_vline(xintercept = 0,linetype=2,color="grey")+
    #geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_shape_manual(values = c(17,16))+
    labs(title = title)+
    geom_segment(data = db_rda_forward_pr.env, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'black') +
    #geom_segment(data = db_rda_forward_pr.spe, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'grey') +
    #geom_text(data = db_rda_forward_pr.spe, aes(CAP1 * 1.1, CAP2 * 1.1, label = name), color = 'grey', size = 3) +
    # geom_text_repel(data=db_rda_forward_pr.spe,aes(x = CAP1*1.1, 
    #                                                y = CAP2 * 1.1, 
    #                                                label = name),size=3,
    #                 color = "grey") +
    geom_text_repel(data=db_rda_forward_pr.env,aes(x = CAP1, 
                                                   y = CAP2, 
                                                   label = name),
                    size=5,,min.segment.length = unit(0, 'lines'),
                    color = "black") +
    #theme_bw() +
    theme_prism(base_size = 16,border = T)+
    annotate("text", x = coord_x, y = coord_y, label= do.call(expression, nonmetric_label),hjust = 1, vjust = 1,
             size = 16/.pt) +
    #geom_text(data = db_rda_forward_pr.env, aes(CAP1 * 1.1, CAP2 * 1.1, label = name), color = 'blue', size = 3) +
    stat_ellipse(aes(group = Gradient,color=Gradient),
                 level=0.68,alpha=1,
                 type = "norm")+
    scale_color_manual(name = "Gradient",values = envPalette)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  p2
  # Unconstrained variation showed in PCoA ----------
  res_forward_pr.site <- data.frame(db_rda_forward_pr$CA$u)[1:2]
  res_forward_pr.site$name <- rownames(res_forward_pr.site)
  res_forward_pr.site$Site <- pcoa$Site
  res_forward_pr.site$Season <- pcoa$Season
  res_forward_pr.site$Gradient <- pcoa$Gradient
  rda1_exp <- paste('PCOA1:', round((db_rda$CA$eig/sum(db_rda$CA$eig))[1]*100, 2), '%')
  rda2_exp <- paste('PCOA2:', round((db_rda$CA$eig/sum(db_rda$CA$eig))[2]*100, 2), '%')
  title = sprintf("Unconstrained ordination (%s) (%s)",name,"bray curtis distance")
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
  
  p3
  # Calculate environment parameter contribution. -------------
  levels <- db_rda_forward_pr$anova %>% 
    as.data.frame() %>% 
    rownames() %>%
    str_remove("\\+ ") %>%
    str_remove_all('`') %>%
    head(-1)
  Contribut_df <- merge_tab_1 %>%
    relocate(levels, .before = 1)
  R_square_table <- data.frame()
  #rm(formula_cum)
  for (i in colnames(Contribut_df)) {
    i1 <- i
    i <- paste0("`",i,"`")
    formula <- as.formula(paste0("taxon~",i))
    set.seed(123)
    uni_dbrda <- capscale(formula, merge_tab_2, distance = 'bray')
    uni_i <- max(RsquareAdj(uni_dbrda)$adj.r.squared,0)
    uni_p <- anova(uni_dbrda)$`Pr(>F)`[1]
    if (!exists("formula_cum")) {
      formula_cum <- formula
    }else{
      formula_cum <- as.formula(paste0(Reduce(paste, deparse(formula_cum)),"+",i))
    }
    set.seed(123)
    db_rda <- capscale(formula_cum,
                       data = merge_tab_2,dist = "bray")
    cum_i <- RsquareAdj(db_rda)$adj.r.squared
    re_tb <- data.frame(name = i1, uni = uni_i,unip  = uni_p, cum = cum_i)
    R_square_table <- rbind(R_square_table,re_tb)
  }
  R_sign_table <- R_square_table %>% 
    mutate(BHP = p.adjust(unip,method = "BH")) %>% 
    filter(BHP < 0.05)
  
  return_result <- list(pcoa=p1,dbRDA=p2,pcoa_dbRDA=p3,R_sign_table = R_sign_table,
                        cumRdf = R_square_table,levels = levels)
  return(return_result)
}
### Function Barplot -----------------------------------------------------------
Barplot <- function(data,sign_v)
{
  # data <- Species_list$cumRd
  # sign_v <- Species_list$levels
  data <- data %>%
    mutate(name = factor(name,levels = rev(name))) %>%
    pivot_longer(2:4,names_to = "group",values_to = "value")
  sign_R <- data %>%
    filter(name %in% tail(sign_v,1),group %in% "cum") %>%
    .$value
  data <- data %>%
    filter(!group == "unip")
  p3 <- ggplot(data) +
    geom_bar(stat = "identity",position = "dodge",aes(x=name,y=value,fill = group,group = group)) +
    scale_fill_manual(values = c("grey","black"),labels = c(bquote("Cumulative adjust" ~ R^2),bquote("Univariate adjust" ~ R^2))) +
    coord_flip() + geom_hline(yintercept = sign_R,linetype=2,color="darkred")+
    ylab("Effect size") + xlab("") +
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
Animal_species <- Animal_OTU_df %>%
  #filter(SampleID %in% rownames(COG_path)) %>%
  column_to_rownames("SampleID") %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  left_join(Animal_TAXA,by = c("rowname"= "ID")) %>%
  filter(!is.na(Species)) %>%
  group_by(Species) %>%
  summarise_if(is.double,max) %>%
  column_to_rownames("Species") %>%
  t()  %>%
  as.data.frame() %>%
  rownames_to_column("ID")
write_tsv(Animal_species,"~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/species.tsv")
#Superkingdom_list <- Beta_diversity(paste0(data_path,"superkingdom.otu"),"Superkingdom")
# Phylum_list <- Beta_diversity(paste0(data_path,"phylum.tsv"),"Phylum")
# Class_list <- Beta_diversity(paste0(data_path,"class.tsv"),"Class")
# Order_list <- Beta_diversity(paste0(data_path,"order.tsv"),"Order")
# Family_list <- Beta_diversity(paste0(data_path,"family.tsv"),"Family")
# Genus_list <- Beta_diversity(paste0(data_path,"genus.tsv"),"Genus")
rm(formula_cum)
Species_list <- Beta_diversity(Animal_species,"Species")
#### Species level -------------------------------------------------------------
bar_plot <- Barplot(Species_list$cumRdf,Species_list$levels)
Species_list$dbRDA_withbarplot <- Species_list$dbRDA + bar_plot +
  plot_layout(design = "AAAB") +
  plot_layout(guides = 'collect')  &
  theme(legend.position='bottom') 
pdf("~/Desktop/marine_sediment/scripts_manuscript_new/shelby_data/dbRDA_small_animal_species.pdf",width = 18,height = 10)
Species_list$dbRDA_withbarplot
dev.off()
### retuen part ----------------------------------------------------------------
return_list <- list(#Superkingdom = Superkingdom_list,
                    Phylum = Phylum_list,
                    Class = Class_list,
                    Order = Order_list,
                    Family = Family_list,
                    Genus = Genus_list,
                    Species = Species_list
                    )
return(return_list)
}
#beta_diversity_dbRDA_list <- beta_diversity_dbRDA()
#save(beta_diversity_dbRDA_list, file="~/Desktop/marine_sediment/RDS_files/beta_diversity_dbRDA_2019morevariables_list.RData")
# library(patchwork)
# pdf("~/Desktop/marine_sediment/results/dbRDA_taxon_2019/dbRDA_taxon_2019.pdf",
#     width = 10,height = 30)
#   beta_diversity_dbRDA_list$Phylum$dbRDA_withbarplot + beta_diversity_dbRDA_list$Class$dbRDA_withbarplot #+
#     #beta_diversity_dbRDA_list$Order$dbRDA_withbarplot + beta_diversity_dbRDA_list$Family$dbRDA_withbarplot + beta_diversity_dbRDA_list$Genus$dbRDA_withbarplot + beta_diversity_dbRDA_list$Species$dbRDA_withbarplot + plot_layout(ncol = 1) 
# dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_phylum.pdf",width = 16,height = 10)
Phylum_list$dbRDA_withbarplot
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_class.pdf",width = 16,height = 10)
Class_list$dbRDA_withbarplot
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_order.pdf",width = 16,height = 10)
Order_list$dbRDA_withbarplot
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_family.pdf",width = 16,height = 10)
Family_list$dbRDA_withbarplot
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_genus.pdf",width = 16,height = 10)
Genus_list$dbRDA_withbarplot
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/dbRDA_taxon_2019_species.pdf",width = 16,height = 10)
Species_list$dbRDA_withbarplot
dev.off()

pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_phylum.pdf",width = 10,height = 10)
Phylum_list$pcoa
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_class.pdf",width = 10,height = 10)
Class_list$pcoa
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_order.pdf",width = 10,height = 10)
Order_list$pcoa
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_family.pdf",width = 10,height = 10)
Family_list$pcoa
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_genus.pdf",width = 10,height = 10)
Genus_list$pcoa
dev.off()
pdf("~/Desktop/marine_sediment/results/Shelby_data/PcoA_2019_species.pdf",width = 10,height = 10)
Species_list$pcoa
dev.off()