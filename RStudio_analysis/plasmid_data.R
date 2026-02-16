library(tidyverse)
#please note that the plasmid.csv file was produced in Excel manually, combining both annotations 
#from prokka and BLAST 
plasmid <- read.csv("plasmid.csv")
plasmid <- plasmid %>% mutate(sample=recode(sample,
                                            "56" = "1h RNApolyAdep",
                                            "57" = "24h RNApolyAdep", 
                                            "528" = "1h RNAdep",
                                            "529" = "24h RNAdep"))
plasmid <- plasmid %>% mutate(improved_annotation=recode(name,
                                          "PROKKA_1" = "putative pGP5-D",
                                          "PROKKA_2" = "putative CT583",
                                          "PROKKA_3" = "putative pGP7-D",
                                          "PROKKA_4" = "tyrosine recombinase XerC",
                                          "PROKKA_5" = "COG0305 DNA helicase",
                                          "PROKKA_6" = "putative DUF5597-containing protein",
                                          "PROKKA_7" = "pGP3-D",
                                          "PROKKA_8" = "putative pGP4-D"))
plasmid_table <- plasmid %>% select(name, improved_annotation, sample, mapped) %>% 
  pivot_wider(names_from=sample, values_from=mapped)
colnames(plasmid_table)[1:2] <-c("PROKKA_ID", "Improved_Annotation")
plasmid_table

write.csv(plasmid_table, "plasmid_counts_table.csv", row.names = FALSE)

