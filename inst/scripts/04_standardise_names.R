library("tidyverse")
library("AbNames")


citeseq <- read_delim("../inst/extdata/ADT_clones/merged_adt_clones.csv")
std <- matchToCiteseq(citeseq)

# Standard names use first antigen per study, replace missing rows and fill
citeseq <- dplyr::bind_rows(std, dplyr::anti_join(citeseq, std)) %>%
    dplyr::group_by(Antigen, Study) %>%
    dplyr::fill(Antigen_std)


# antijoin citeseq and std
#



# Check for entries with no data name match

#citeseq %>%
#    filter(is.na(ALT_ID), ! Isotype_Control) %>%
#    select(Antigen, Data_Name, Study, Cat_Number, Clone)


