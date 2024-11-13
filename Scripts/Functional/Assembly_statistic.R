library(tidyverse)
Assembly_df <- read_csv("~/Desktop/marine_sediment/manuscript/Figure_v5/Assembly_table.csv")
N50_df <- read_tsv("~/Desktop/marine_sediment/manuscript/Figure_v5/N50.txt")
df <- Assembly_df %>%
  select(Name,`Contig number`,`Genome size (bp)`) %>%
  left_join(N50_df)

library(xlsx)
write.xlsx(df, "~/Desktop/marine_sediment/manuscript/manuscripts/Marine_microbiom_resubmit_Cell_Reports/Supplementary_tables_new/Supplementary_Tabs.xlsx",
           sheetName="Table S5",row.names = T,append=TRUE)

#write_tsv(df, "~/Desktop/marine_sediment/manuscript/Figure_v5/tables/TabS7_assembly_summary.tsv")
