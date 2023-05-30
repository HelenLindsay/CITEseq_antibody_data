# This script matches antibody names with identifiers using the gene_aliases
# table in the package AbNames

library("tidyverse")
library("AbNames")
library("janitor")
library("stringi")

# To do: group and check differences in ALT_ID.
# Chung c-Fos has wrong catalogue number?

existing <- ls()

if(! grepl("vignettes", getwd())) { setwd("./vignettes") }

#merged_clones <- "../inst/extdata/ADT_clones/merged_adt_clones.tsv"
#citeseq <- readr::read_delim(file = merged_clones)

# Make a table for querying alias table, excluding isotype controls ----

# Remove non-ascii characters
citeseq <- all_clones %>%
    dplyr::mutate(across(where(is.character),
                         ~stringi::stri_trans_general(.x,
                                id="Any-Latin;Greek-Latin;Latin-ASCII")))

# Add an ID by pasting study/antigen
citeseq <- AbNames::addID(citeseq)
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen, Isotype_Control)

query_df <- AbNames::makeQueryTable(citeseq_q,
                                    ab = "Antigen",
                                    control_col = "Isotype_Control")

# To allow re-running, remove ID columns before re-annotating -----
id_cols <- c("ALT_ID", "HGNC_ID", "HGNC_SYMBOL", "ENSEMBL_ID",
             "UNIPROT_ID", "ENTREZ_ID", "SOURCE")

citeseq <- all_clones %>%
    dplyr::select(-any_of(id_cols))

# Add an ID by pasting study/antigen
citeseq <- AbNames::addID(citeseq)
citeseq_q <- citeseq %>% dplyr::select(ID, Antigen, Isotype_Control)

query_df <- AbNames::makeQueryTable(citeseq_q,
                                    ab = "Antigen",
                                    control_col = "Isotype_Control")

# Annotate using the gene aliases table ----
alias_results <- searchAliases(query_df)

# Remove matches to several genes, select just columns of interest
alias_results <- alias_results %>%
    # Select ID and HGNC columns
    dplyr::select(matches("ID|HGNC"), name, SOURCE) %>%
    unique() %>% # Collapse results with same ID from different queries
    group_by(ID) %>%
    # Collapse multi-subunit entries, convert "NA" to NA
    dplyr::summarise(dplyr::across(all_of(id_cols), ~toString(unique(.x)))) %>%
    dplyr::mutate(dplyr::across(all_of(id_cols), ~na_if(.x, "NA")))

# Re-add entries into citeseq table
citeseq <- citeseq %>%
    dplyr::left_join(alias_results, by = "ID") %>%
    dplyr::relocate(ID, Study, Antigen, Cat_Number, HGNC_ID, ALT_ID) %>%
    unique()

# Annotate with BioLegend information using catalogue number or clone
citeseq <- searchTotalseq(citeseq)

# Check for missing annotation ----
citeseq %>%
    dplyr::filter(is.na(ALT_ID) & ! Isotype_Control)

# Patch the antigens that are still missing -----
cs_patch <- tibble::tribble(~Antigen, ~value,
                            "Annexin V", "annexin V",
                            "DopamineD4receptor", "dopamine receptor D4",
                            "DopamineReceptorD4", "dopamine receptor D4",
                            # There are antibodies against CD235a and CD235a/b
                            "CD235", "CD235a",
                            "CD32/ Fcg RII", "CD32",
                            "CD110 MPL", "CD110")

cs_patch <- cs_patch %>%
    dplyr::left_join(gene_aliases, by = "value") %>%
    dplyr::select(any_of(colnames(citeseq)))

citeseq <- citeseq %>%
    dplyr::rows_patch(cs_patch) %>%
    dplyr::ungroup()

# Fill Cat_Number, Oligo_ID, etc ----

# First fill other entries given (TotalSeq) catalogue number:
# If catalogue number is available, can fill clone, oligo
# For BD Bioscieces, I'm not sure if Antigen and Clone is enough to fill
# Cat_Number

exp_nrow <- nrow(citeseq)
n_na <- citeseq %>%
    dplyr::summarise(across(c("Cat_Number", "Clone", "Oligo_ID"),
                            ~sum(is.na(.x))))

citeseq_custom <- citeseq %>%
    dplyr::filter(grepl("custom", Cat_Number) | isTRUE(Custom_Antibody))


# TO DO: CHECK FOR ERRORS AND REMOVE "IGNORE"
citeseq <- citeseq %>%
    dplyr::filter(! grepl("custom", Cat_Number) &
                      (is.na(Custom_Antibody) | ! Custom_Antibody)) %>%
    fillByGroup(group = "Cat_Number",
                fill = c("Clone", "Oligo_ID", "TotalSeq_Cat", "ALT_ID",
                         "HGNC_ID", "ENSEMBL_ID"),
                multiple = "ignore") %>%
    # Allow leniency for Reactivity as it hasn't been sorted
    fillByGroup(group = "Cat_Number", fill = c("Reactivity"),
                multiple = "ignore") %>%
    fillByGroup(group = "RRID", fill = c("Clone", "Reactivity"),
                multiple = "ignore") %>%
    # It is possible to share oligo, antigen and clone but not cat number
    fillByGroup(group = c("Clone", "Oligo_ID", "TotalSeq_Cat"),
                fill = c("Cat_Number"), multiple = "ignore")

citeseq <- dplyr::bind_rows(citeseq, citeseq_custom) %>%
    dplyr::arrange(Study, Antigen)

exp_nrow == nrow(citeseq)

print(n_na)

citeseq %>%
    dplyr::summarise(across(c("Cat_Number", "Clone", "Oligo_ID"),
                            ~sum(is.na(.x))))


# Save annotated data set ----
merged_clones <- "../inst/extdata/ADT_clones/merged_adt_clones.tsv"
readr::write_delim(citeseq, file = merged_clones)

rm(list = setdiff(c(existing, "citeseq")))

# Hao_2021 and Liu_2021 have entries for some but not all Cat_Numbers
# It appears info can be added for Liu, but Hao entries are custom
# (But do not ever have other group members)


# Standardising names for controls
#dplyr::group_by(Cat_Number) %>%
#dplyr::mutate(n = ifelse(n == 1 & grepl("kappa", Suggested_Antigen),
#                         100, n)) %>% # Prefer kappa if it's there

# HGNC:12102 = TRAV1-2 = TCR Va7.2?
# TRAV24, TRAJ18 = TCRVa24-Ja18
# HLA-DR = HLA-DRA?

# For antigens like KIR2DL1/S1/S3/S5, want to split like subunit
# CD18 associates with CD11a (LFA-1) CD11-b (Mac-1)
# CD3 (CD3E)__Nathan_2021 should match CD3D,E and G?
