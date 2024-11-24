# Title: Quality control
# Project: Chemical Pollution Drives Taxonomic and Functional Shifts in Marine Sediment Microbiome, Influencing Benthic Metazoans
# Maintainer: Lin-Lin Xu
# Last update: Nov 24, 2024
# License: MIT

setwd("~/Desktop/marine_sediment/results/")
library(ggplot2)
library(reshape2)
library(tidyverse)
Final_ReadSummaryStatistics <- read.delim("~/Desktop/marine_sediment/results/Final_ReadSummaryStatistics.txt")
test <- melt(Final_ReadSummaryStatistics,id=c("sample","human_count","Adapter","Duplication","perc_usable"))
test1 <- melt(Final_ReadSummaryStatistics,id=c("sample","raw_PE_count","usable_read"))
## Huamn percentage 
df_percentage <- Final_ReadSummaryStatistics 
library(xlsx)
write.xlsx(df_percentage, "~/Desktop/marine_sediment/manuscript/manuscripts/Supplementary_tables_new/Supplementary_Tabs.xlsx",
           sheetName="Table S2",row.names = F,append=TRUE)
###bar plot###
ggplot(test,aes(x=sample,y=value,fill=variable))+geom_bar(stat='identity', position='dodge') + theme_bw() + 
  labs(x = 'Sample',y = 'Reads counts') + theme(axis.text.x = element_text(angle = 90,hjust = 1))
###stack plt###
QC_P <- ggplot(test1, aes(sample,value, fill=variable))+
  geom_bar(stat='identity',position='fill') + labs(x = '',y = 'Percentage') + theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) + labs(fill="")
###circle plot###
ggplot(test1, aes(sample,value, fill=variable))+
  geom_bar(stat='identity',width=0.5,position='stack',size = 5) +
  coord_polar("y", start=0)
  #geom_text(aes(y=value,label=value),colour="white")
result_list <- list(df = Final_ReadSummaryStatistics,
                    P = QC_P)
saveRDS(result_list,"~/Desktop/marine_sediment/RDS_files/QC.rds")
