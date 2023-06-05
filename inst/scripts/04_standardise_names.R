library("tidyverse")
library("AbNames")




citeseq <- read_delim("inst/extdata/ADT_clones/merged_adt_clones.tsv")


# Check for entries with no data name match

#citeseq %>%
#    filter(is.na(ALT_ID), ! Isotype_Control) %>%
#    select(Antigen, Data_Name, Study, Cat_Number, Clone)


