# Notes ----
# This script formats the antibody clone tables from each data set
# into a single table wtih consistently named columns

# Common column names
# Antigen
# Clone
# Cat_Number: catalogue number
# Study: firstauthor_year (we give more descriptive names later)
# RRID: Research Resource ID
# TotalSeq_Cat: A/B/C for different sequencing protocols
# Oligo_ID: Just the numeric part of the tag e.g. 0034 for A0034
# Ensembl_ID
# Reactivity: Which species does the antibody react against? Usually human
# Origin_Species: Which species was the antibody generated in?
# Vendor: Supplier of the antibody
# Dilution: how to write this in a standard way? In arunachalam, buus
# Barcode_Sequence: sequence of the oligo
# Isotype: the isotype of the antibody
# Isotype_Control (logical): is the antibody an isotype control?

# Additional columns are specific to one study.
# Column "Gene_Name" is used when a gene symbol is provided.
# At this stage in the data preparation we do not check if this is the
# official gene symbol.

# Code ----

existing <- ls()
if(! grepl("vignettes", getwd())) { setwd("./vignettes") }

library("dplyr")
library("readr")
library("tidyr")
library("stringr")
library("AbNames")

clone_dir <- "../inst/extdata/ADT_clones"

# Read in the ADT clone tables ----
# TotalC is a table of a standard antibody panel from BioLegend

adt_clone_fnames <- list.files(clone_dir, full.names = TRUE)

adt_clone_fnames <- adt_clone_fnames[! grepl("TotalSeq_C|merged",
                                             adt_clone_fnames)]
adt_clones <- lapply(adt_clone_fnames, readr::read_delim)
names(adt_clones) <- gsub("\\..*", "", basename(adt_clone_fnames))

totalseqC <- list.files(clone_dir, full.names = TRUE, pattern ="^TotalSeq_C")

tenx_clone_fnames <- list.files("../inst/extdata/10x_antibody_info",
                                full.names=TRUE)
tenx_clones <- lapply(tenx_clone_fnames, readr::read_delim)
names(tenx_clones) <- gsub("\\..*", "", basename(tenx_clone_fnames))

tenx_id_to_fname <- tenx_clones$`10x_files`
tenx_clones$`10x_files` <- NULL


# Arunachalam_2020 ----
arunachalam <- adt_clones$Arunachalam_2020
arunachalam <- arunachalam %>%
    dplyr::rename(Cat_Number = `Cat #`,
                  Vendor = Supplier,
                  TotalSeq_Tag = Tag) %>%
    dplyr::mutate(Study = "Arunachalam_2020",
                  Cat_Number = as.character(Cat_Number),
                  Vendor = ifelse(Vendor == "Biolegend",
                                  "BioLegend", Vendor),
                  Dilution = gsub("^0|:00$", "", as.character(Dilution))) %>%
    dplyr::select(- `...7`) %>%
    tidyr::separate(TotalSeq_Tag, into = c("TotalSeq_Cat", "Oligo_ID"),
                    sep = "(?<=[ABC])") %>%
    dplyr::mutate_if(is.factor, as.character)

# Bai_2022 ----

bai <- adt_clones$Bai_2022 %>%
    dplyr::mutate(Study = "Bai_2022")

# Buus_2021 ----

buus <- adt_clones$Buus_2021
buus <- buus %>%
    dplyr::rename_with(function(x) gsub("\n", "", x)) %>%
    dplyr::rename(Antigen= Marker,
                  Isotype = `Isotype(Mouse)`,
                  Cat_Number = `Cat#`,
                  Oligo_ID = `TotalSeq CTag`) %>%
    dplyr::mutate(Study = "Buus_2021",
                  Origin_Species = ifelse(is.na(Isotype), NA, "Mouse"),
                  TotalSeq_Cat = "C",
                  Oligo_ID =
                      stringr::str_pad(as.character(Oligo_ID), 4, pad = "0"),
                  Cat_Number = as.character(Cat_Number)) %>%
    # Drop concentrations and Alias (not needed to identify antibody)
    dplyr::select(-matches("conc"), -Alias) %>%
    dplyr::mutate_if(is.factor, as.character)


# Cadot_2020  ----

# TotalSeq_Cat added by looking up catalog numbers
# Catalogue number appears incorrect for CD49
cadot <- adt_clones$Cadot_2020 %>%
    dplyr::mutate(Study = "Cadot_2020",
                  Vendor = "BioLegend",
                  TotalSeq_Cat = "A",
                  Catalog = as.character(Catalog)) %>%
    dplyr::rename(Antigen = Antibody,
                  Cat_Number = Catalog) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::select(-Stain)

# Chung_2021 ----
# Method inCITE-seq = intranuclear protein and transcriptome
# Totalseq A format

chung_2021 <- tibble::tribble(~Antigen, ~Clone, ~Cat_Number,
                              "p65", "Poly6226", "622601")#,
# Used in mouse hippocampus
#"NeuN", "1B7", "834502",
#"PU.1", "7C2C34", "681307",
#"c-Fos", "Poly6414", "834502")

chung_2021 <- chung_2021 %>%
    dplyr::mutate(Study = "Chung_2021",
                  Vendor = "BioLegend",
                  TotalSeq_Cat = "A")

# Colpitts_2023 ----
colpitts <- adt_clones$Colpitts_2023 %>%
    # Make concentration consistent with other studies
    dplyr::mutate(ug_per_100_uL = 10*ug_per_100_uL) %>%
    dplyr::rename(Conc_ug_per_ml = ug_per_100_uL)

# Fernandez_2019 ----
fernandez <- adt_clones$Fernandez_2020
fernandez <- fernandez %>%
    dplyr::rename(Antigen = Name,
                  Vendor = Company,
                  Cat_Number = `Catalog Number`) %>%
    dplyr::mutate(Study = "Fernandez_2019",
                  Cat_Number = as.character(Cat_Number),
                  Vendor = ifelse(Vendor == "Biolegend",
                                  "BioLegend", Vendor)) %>%
    dplyr::select(Antigen, Clone, Vendor, Cat_Number, Study) %>%
    dplyr::mutate_if(is.factor, as.character)

# Fidanza 2020 ----
fidanza <- adt_clones$Fidanza_2020 %>%
    dplyr::mutate(TotalSeq_Cat = gsub("TotalSeq-", "", Antibody_Cat)) %>%
    dplyr::select(-Antibody_Cat) %>%
    dplyr::rename(Conc_ug_per_ml = `Concentration(migrograms/ml)`)

# Frangieh 2021 ----

frangieh <- adt_clones$Frangieh_2021 %>%
    dplyr::rename(Antigen = `Target `,
                  Clone = `Clone `,
                  Barcode_Sequence = Barcode) %>%
    dplyr::mutate(Oligo_ID = sprintf("%04d", `TotalSeq A#`),
                  TotalSeq_Cat = "A",
                  Vendor = "BioLegend",
                  Study = "Frangieh_2021") %>%
    dplyr::select(Antigen, Clone, Barcode_Sequence, Oligo_ID,
                  Vendor, Study)


# Granja_2019 ----

granja <- adt_clones$Granja_2019
granja <- granja %>%
    dplyr::rename(Antibody_Category = Category,
                  Cat_Number = Catalog,
                  Antigen = Specificity,
                  TEMP = Name, # This might be the name used in the data
                  TotalSeq_Cat = Category,
                  Barcode_Sequence = `Barcode Sequence`,
                  Oligo_ID = Barcode) %>%
    dplyr::mutate(Study = "Granja_2019",
                  Cat_Number = as.character(Cat_Number),
                  TotalSeq_Cat = gsub(".*([ABC])$", "\\1", TotalSeq_Cat),
                  Vendor = "BioLegend") %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::select(-Condition)


# Hao_2021 ----

# Two assays - CITE-seq and ECCITE-seq.  Siglec-8 and CD26 are in both.
# In ECCITE-seq, CD45RA and CD45RO each have two barcodes - a custom preparation
# and TotalSeq C.

# Using Specificity as antigen name instead of Antigen
# Antigen in original table isn't needed for matching data

# Patch clones that got converted to numeric at some point
hao_patch <- tibble(Cat_Number =  c("310729", "328135"),
                    Clone = c("5E8", "5E10"))

hao <- adt_clones$Hao_2021 %>%
    #dplyr::rename(Data_Name = Antigen) %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  Antigen = ifelse(is.na(Specificity),
                                   Antigen, Specificity),
                  # We assume if Clone is identical to Antigen, this is an error
                  Clone = ifelse(Antigen == Clone, NA, Clone)) %>%
    dplyr::rows_update(hao_patch, by = "Cat_Number") %>%
    dplyr::select(-Specificity)


# Holmes_2020 ----
holmes_2020 <- adt_clones$Holmes_2020 %>%
    dplyr::select(-`...6`, -`...7`)

# Kaufmann_2021 ----
kaufmann <- adt_clones$Kaufmann_2021
kaufmann <- kaufmann %>%
    dplyr::rename(Antigen = antigen,
                  Vendor = supplier,
                  Cat_Number = `order no`,
                  Lot = lot,
                  Barcode_Sequence = barcode,
                  TEMP = `marker_ name`, # Check if it matches data
                  Clone = clone) %>%
    dplyr::mutate(Study = "Kaufmann_2021",
                  Cat_Number = as.character(Cat_Number),
                  Vendor = ifelse(Vendor == "Biolegend",
                                  "BioLegend", Vendor)) %>%
    dplyr::mutate_if(is.factor, as.character)


# Kotliarov_2020 -----
# Cat_Numbers with * are custom orders similar to
kotliarov <- adt_clones$Kotliarov_2020
kotliarov <- kotliarov %>%
    dplyr::rename(TotalSeq_Cat = Category,
                  Antigen = Specificity,
                  Cat_Number = `Cat. Num.`,
                  Reactivity = Species,
                  Oligo_ID = Barcode_ID,
                  Custom_Antibody = `Custom `,
                  Barcode_Sequence = Barcode) %>%
    dplyr::mutate(Study = "Kotliarov_2020",
                  Cat_Number = as.character(Cat_Number),
                  TotalSeq_Cat = gsub("TotalSeq", "", TotalSeq_Cat),
                  Oligo_ID = stringr::str_pad(as.character(Oligo_ID),
                                              4, pad = "0"),
                  Vendor = "BioLegend") %>%
    dplyr::mutate_if(is.factor, as.character)

# Krebs 2020 ----
# No clone information was given.

# Lawlor_2021 ----
lawlor <- adt_clones$Lawlor_2021
lawlor <- lawlor %>%
    dplyr::rename(Antigen = Antibody,
                  Barcode_Sequence = `CITE-seq barcode`,
                  Vendor = Supplier,
                  Reactivity = Reactive) %>%
    dplyr::select(-Note) %>%
    dplyr::mutate(Study = "Lawlor_2021") %>%
    dplyr::mutate_if(is.factor, as.character)


# Leader_2021 ----

nm <- "Leader_2021"

leader <- adt_clones$Leader_2021 %>%
    dplyr::mutate(Oligo_ID = AbNames:::.noDups(
        gsub("TotalSeq.*-[AC]([0-9]+).*", "\\1", Antigen), Antigen),
        TotalSeq_Cat = AbNames:::.noDups(
            gsub("TotalSeq.*-([AC])[0-9]+.*", "\\1", Antigen), Antigen),
        Antigen_cp = gsub(" .*", "", Antigen),
        Antigen = gsub(".*[Aa]nti-(.*)", "\\1", Antigen),
        Antigen = gsub(" Antibody", "", Antigen),
        Study = "Leader_2021",
        RRID = AbNames:::.noDups(
            gsub(".*RRID:(.*)$", "\\1", Cat_Number), Cat_Number),
        Cat_Number = gsub(".*#[: ]*([^;]+).*", "\\1", Cat_Number),
        Cat_Number = as.character(Cat_Number),
        Clone = AbNames:::.noDups(
            gsub(".*\\(clone (.*)\\)", "\\1", Antigen), Antigen),
        Reactivity = AbNames:::.noDups(
            gsub("([Hh]uman[^ ]*) .*", "\\1", Antigen), Antigen),
        Antigen = gsub("([Hh]uman[^ ]*) (.*)", "\\2", Antigen),
        Antigen = gsub(" \\(clone .*\\)", "", Antigen),
        Antigen = ifelse(Antigen == "purified",
                         Antigen_cp, Antigen),
        # Antigen "CD8 – 146Nd" causes problems, I don't know why
        Antigen = ifelse(grepl("CD8", Antigen), "CD8", Antigen)) %>%
    dplyr::select(-Antigen_cp)

# LeCoz_2021 ----

# Totalseq C universal cocktail, 132 + 8 isotype control
total_c <- read_delim(paste0("../inst/extdata/ADT_clones/",
  "TotalSeq_C_Human_Universal_Cocktail_v1_137_Antibodies_399905_Barcodes.csv"))

lecoz <- total_c %>%
  dplyr::rename(Barcode_Sequence = Barcode) %>%
  dplyr::mutate(Vendor = "BioLegend",
                TotalSeq_Cat = "C",
                Oligo_ID = gsub("C", "", DNA_ID),
                Study = "LeCoz_2021",
                Reactivity =
                  gsub("anti-([a-z\\/]+) .*", "\\1", Description),
                Antigen = gsub("anti-[a-z\\/]+ (.*)", "\\1", Description),
                Reactivity = AbNames:::.noDups(Reactivity, Antigen)) %>%
  dplyr::select(Antigen, Clone, Oligo_ID, Vendor, TotalSeq_Cat, Study,
                Reactivity)

# Lee_2021 ----
lee <- adt_clones$Lee_2021
lee <- lee %>%
    dplyr::mutate(Study = "Lee_2021",
                  Cat_Number = as.character(Cat_Number),
                  TotalSeq_Cat = gsub("^TotalSeq-([ABC])[0-9]+$",
                                      "\\1", Antibody_Category),
                  Oligo_ID = gsub("^.*[ABC]([0-9]+)$", "\\1",
                                  Antibody_Category)) %>%
    dplyr::mutate(Vendor = "BioLegend",
                  Clone = ifelse(grepl("clone", Antigen),
                                 gsub(".*\\(clone: (.*)\\)$", "\\1", Antigen),
                                 NA),
                  Antigen = gsub("\\(clone.*$", "", Antigen)) %>%
    dplyr::select(-Antibody_Category) %>%
    dplyr::mutate_if(is.factor, as.character)

# Li 2023 ----
# clone information manually formatted from Sup Table 1

li_2023 <- adt_clones$Li_2023 %>%
    dplyr::mutate(Dilution = as.character(Dilution),
                  Dilution = gsub("^([0-9]+:[0-9]+):.*", "\\1", Dilution)) %>%
    dplyr::rename(Cat_Number = TotalSeqCat,
                  Antigen = Antibody) %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number))

# Liu_2021 ----

# From reagent table, TotalSeq C custom human panel was used

# Here the TotalSeq suffix is called DNA_ID
liu <- adt_clones$Liu_2021
liu <- liu %>%
    dplyr::mutate(Study = "Liu_2021",
                  Vendor = "BioLegend",
                  Antigen =
        gsub("^[0-9]+ (anti-mouse/human)?(anti-[hH]uman)?(/mouse)?(/rat)? ?",
             "", Description),
                  Reactivity =
        gsub("[0-9]+ ((anti-mouse/human)?(anti-[hH]uman)?(/mouse)?(/rat)?).*",
                                    "\\1", Description),
                  RRID =
        gsub("(AB_[0-9]+).*", "\\1", `RRID (Biolegend Catalogue Number)`),
                  Cat_Number =
        gsub('^.* Cat. No. ([0-9]+)).*$', "\\1",
             `RRID (Biolegend Catalogue Number)`),
             Oligo_ID = gsub("^[ABC]", "", DNA_ID),
             TotalSeq_Cat = gsub("^([ABC])[0-9]+$", "\\1", DNA_ID))%>%
    dplyr::relocate(Antigen) %>%
    dplyr::select(-`RRID (Biolegend Catalogue Number)`,
                  -Description, -DNA_ID) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  Reactivity = gsub("anti-", "", Reactivity)) %>%
    dplyr::rename(Barcode_Sequence = Barcode)

# Mair_2020 ----

# Selected the -Ab-0 from the full clone table

mair <- adt_clones$Mair_2020
mair <- mair %>%
    dplyr::rename(Vendor_Cat_RRID = RRID) %>%
    dplyr::mutate(Antigen = gsub("-Ab-O", "", Antigen),
                  Clone = paste0(Clone, Cat_Number),
                  Clone = gsub("^.*\\(clone(.*)\\).*$", "\\1", Clone),

                  Vendor_Cat_RRID = gsub(";.*$", "", Vendor_Cat_RRID),
                  Oligo_ID =
                      gsub("^.*\\s+(Cat#)?([^ ]+)$", "\\2", Vendor_Cat_RRID),
                  Oligo_ID = ifelse(Oligo_ID == "core", NA, Oligo_ID),

                  Vendor = "BD Biosciences",
                  Study = "Mair_2020") %>%

    dplyr::select(Antigen, Clone, Oligo_ID, Vendor, Study)

# Mimitou_2019 ----
mimitou_2019 <- adt_clones$Mimitou_2019
mimitou_2019 <- mimitou_2019 %>%
    dplyr::rename(Antigen = Antibody,
                  Clone = clone,
                  Barcode_Sequence = barcode) %>%
    dplyr::mutate(Study = "Mimitou_2019")

# Mimitou_2021_cite_and_dogma_seq ----

mimitou_2021 <- adt_clones$Mimitou_2021_cite_and_dogma_seq
mimitou_2021 <- mimitou_2021 %>%
    dplyr::rename(Barcode_Sequence = Barcode) %>%
    dplyr::mutate(TotalSeq_Cat = gsub("^.*-([ABC])[0-9]+$", "\\1", Name),

                  Oligo_ID = gsub("^.*-[ABC]([0-9]+)", "\\1", Name),

                  Reactivity = gsub("^anti-([^ ]+)\\s.*$", "\\1", Specificity),
                  Reactivity =
                     ifelse(grepl("Rat", Specificity), NA, Reactivity),

                  Antigen = gsub("^anti-[^ ]+\\s(.*)$", "\\1", Specificity),

                  DOGMAseq_Panel = is.na(`not present in the DOGMA-seq panel`),

                  Study = "Mimitou_2021",

                  Vendor = ifelse(grepl("TotalSeq", Name),
                                  "BioLegend", "NA")) %>%
    dplyr::select(Antigen, TotalSeq_Cat, Clone, Vendor,
                  Oligo_ID, Barcode_Sequence, Study)
   # TO DO: COLLECT DOGMA-SEQ INTO ASSAY COLUMN TO MATCH HAO

# Nathan_2021 ----

# Useful columns Biological relevance and Titration to put into sample metadata

nathan <- adt_clones$Nathan_2021
nathan <- nathan %>%
    dplyr::rename(Antigen = Marker,
                  Clone = `Clone ID`,
                  Cat_Number = `Catalogue #`) %>%
    dplyr::select(Antigen, Clone, Cat_Number) %>%
    dplyr::mutate(Study = "Nathan_2021",
                  Cat_Number = as.character(Cat_Number))

# Nettersheim_2022 ----

# Preprocessed from supplementary table S1 by matching to TotalSeq C
# universal set
nettersheim <- adt_clones$Nettersheim_2022 %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  # Thr181 refers to phospho Tau (Thr181)
                  Antigen = gsub("Thr181", "Tau Phospho (Thr181)", Antigen)) %>%
    dplyr::select(-Universal)

# Papalexi_2021 -----

papalexi <- adt_clones$Papalexi_2021 %>%
    dplyr::mutate(Study = "Papalexi_2021",
                  Oligo_ID = stringr::str_pad(as.character(Oligo_ID),
                                              4, pad = "0"),
                  across(where(is.character), ~stringr::str_squish(.x))) %>%
    dplyr::rename(Assay = Experiment)

# Pei_2020 ----

pei_2020 <- adt_clones$Pei_2020
pei_2020 <- pei_2020 %>%
    tidyr::separate(Product, sep = " ", extra = "merge",
                    into = c("TotalSeq_Cat", "Oligo_ID",
                             "Reactivity", "Antigen_Long")) %>%
    dplyr::mutate(TotalSeq_Cat = gsub("TotalSeqTM-", "", TotalSeq_Cat),
                  Reactivity = gsub("anti-", "", Reactivity),
                  Study = "Pei_2020",
                  Cat_Number = as.character(Cat_Number)) %>%
    dplyr::select(-Antigen_Long)

# Poch_2021 ----
poch_2021 <- adt_clones$Poch_2021 %>%
    dplyr::mutate(Study = "Poch_2021",
                  Cat_Number = as.character(Cat_Number))


# PomboAntunes_2021 ----

# Note that when splitting into subexperiments, CD207_mh occurs in both panels

# This includes a human and a mouse panel

# Names for the same isotype control used in different concentrations differ
pombo <- adt_clones$PomboAntunes_2021 %>%
    dplyr::rename(RRID = `Identifier RRID`,
                  Cat_Number = `Cat. No. BioLegend`,
                  Barcode_Sequence = Barcode,
                  Antigen = Name,
                  Gene_Name = Gene) %>%
    tidyr::separate(`TotalSeq ID`, into = c("TotalSeq_Cat","Oligo_ID"),
                    sep = "(?<=[ABC])") %>%
    dplyr::mutate(Study = "PomboAntunes_2021",
                  Vendor = "BioLegend",
                  # Dilution is more commonly formatted as 1:10 than 1/10
                  Dilution = gsub("\\/", ":", Dilution),
                  # Standardise the names of isotype controls
                  Antigen = ifelse(grepl("^IgG.*-", Antigen),
                                   gsub("-", "_", Antigen), Antigen))

# Pont 2020 ----

# Column microgram_per_test_1x could go into sample metadata

pont <- adt_clones$Pont_2020
pont <- pont %>%
    dplyr::rename(Antigen = Marker,
                  Cat_Number = Ref,
                  TotalSeq_Cat = Category,
                  Oligo_ID = Barcode,
                  Barcode_Sequence = BarcodeSequence) %>%
    dplyr::select(-microgram_per_test_1x) %>%
    dplyr::mutate(Study = "Pont_2020",
                  Cat_Number = as.character(Cat_Number),
                  TotalSeq_Cat = gsub("TotalSeq-", "", TotalSeq_Cat))

# Qi 2020 ----

qi_2020 <- adt_clones$Qi_2020
qi_2020 <- qi_2020 %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number))


# RinconArevalo_2021 ----

rincon <- adt_clones$RinconArevalo_2021 %>%
    dplyr::select(-`...7`) %>%
    # Dilution has been imported as a time, fix
    dplyr::mutate(Dilution = gsub(":00$", "", as.character(Dilution)))

# Shangguan_2021 ----

shangguan <- adt_clones$Shanguan_2021
shangguan <- shangguan %>%
    dplyr::rename(Vendor = `Source or reference`,
                  Dilution = `Additional information`) %>%
    dplyr::select(Designation, Identifiers, Vendor, Dilution) %>%
    tidyr::separate(Identifiers, c("Cat_Number","Clone","RRID"),
                    sep = ",[ ]?", fill = "right") %>%
    dplyr::mutate(Origin_Species =
                      gsub("^.*\\(([^ ]+)\\s.*$", "\\1", Designation),
                  Reactivity = case_when(
                      grepl("Anti-Human Mouse", Designation) ~ "mouse",
                      grepl("Anti-Human Rat", Designation) ~ "rat",
                      TRUE ~ "human"),
                  Cat_Number = gsub("Cat# ", "", Cat_Number),
                  Clone = gsub("Clone ", "", Clone),
                  RRID = gsub("RRID:", "", RRID),
                  Antigen = gsub("^Anti-[^ ]+\\s(.*)\\(.*$", "\\1",
                                 Designation),
                  Antigen = stringr::str_squish(Antigen),
                  Study = "Shangguan_2021",
                  TotalSeq_Cat = "C",
                  Dilution = gsub("^([^:]*):[ \n]", "", Dilution)) %>%
    dplyr::select(-Designation)

# Stephenson_2021 ----

# From methods in paper:
# (Newcastle, Cambridge)
# cells were stained with the custom panel TotalSeq-C (BioLegend, 99813)
# (London)
# Cells were stained with TotalSeq-C antibodies (BioLegend, 99814)

# Searching BioLegend website for 99813 or 99814 returns same panel
# Not all antibodies in TotalSeq panel are in Stephenson and vice versa

totalseqC <- readr::read_delim(totalseqC) %>%
    dplyr::rename(Antibody = Description)

stephenson_patch <- data.frame(Clone = "M7004D06",
                               Antigen = "Tau Phospho (Thr181)",
                               Reactivity = "human")

stephenson <- adt_clones$Stephenson_2021
stephenson <- stephenson %>%
    AbNames:::left_join_any(totalseqC,
                            cols = c("Antibody", "Barcode", "Clone")) %>%
    dplyr::rename(Barcode_Sequence = Barcode,
                  Ensembl_ID = `Ensemble ID`,
                  Gene_Name = `Gene name`) %>%
    dplyr::mutate(Reactivity = gsub("^anti-([^ ]+)\\s.*$", "\\1", Antibody),
                  Reactivity = ifelse(Reactivity == Antibody, NA, Reactivity),
                  Antigen =  gsub("^anti-[^ ]+\\s(.*)$", "\\1", Antibody),
                  Study = "Stephenson_2021",
                  Vendor = "BioLegend") %>%
    tidyr::separate(DNA_ID, into = c("TotalSeq_Cat","Oligo_ID"),
                    sep = "(?<=[ABC])") %>%
    dplyr::select(-Antibody) %>%
    dplyr::rows_update(stephenson_patch)

# Stoeckius_2017 ----
stoeckius_2017 <- adt_clones$Stoeckius_2017 %>%
    dplyr::select(Antigen, Clone, Vendor, Reactivity,
                  Barcode_Sequence, Notes) %>%
    dplyr::mutate(Study = "Stoeckius_2017")

# Stoeckius_2018 ----

stoeckius <- adt_clones$Stoeckius_2018
stoeckius <- stoeckius %>%
    dplyr::rename(Antigen = Protein,
                  Vendor = Supplier) %>%
    dplyr::mutate(Study = "Stoeckius_2018")


# Stuart_2019 ----

stuart <- adt_clones$Stuart_2019
stuart <- stuart %>%
    dplyr::rename(Antigen = Name,
                  Cat_Number = `Catalogue #`,
                  Barcode_Sequence = `Barcode Sequence`) %>%
    dplyr::mutate(Study = "Stuart_2019")


# Su_2020 ----

su <- adt_clones$Su_2020
su <- su %>%
    dplyr::rename(Barcode_Sequence = Barcode) %>%
    dplyr::mutate(Description = gsub("^[0-9]*\\s", "", Description),
                  Reactivity = gsub("^anti-([^ ]+)\\s.*", "\\1", Description),
                  Reactivity = ifelse(grepl("Tau|c-Met", Reactivity),
                                      NA_character_, Reactivity),
                  Reactivity = ifelse(Reactivity == Description,
                                      NA_character_, Reactivity),
                  TEMP = ifelse(grepl("sotype", Description), Description,
                                NA_character_),
                  Study = "Su_2020",
                  Antigen = gsub("[hH]uman|anti-|mouse|rat", "", Description),
                  Antigen = gsub("^\\/* ", "", Antigen)) %>%
    tidyr::separate(DNA_ID, into = c("TotalSeq_Cat","Oligo_ID"),
                    sep = "(?<=[ABC])") %>%
    dplyr::mutate(Antigen = ifelse(is.na(TEMP), Antigen, TEMP)) %>%
    dplyr::select(-TEMP)

# Triana 2021 ----

triana_2021 <- adt_clones$Trzupek_2021 %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  Study = "Triana_2021")

# Trzupek_2020 ----

# Trzupek has repeated entries for each subexperiment
trzupek_2020 <- adt_clones$Trzupek_2020
trzupek_2020 <- trzupek_2020 %>%
    dplyr::mutate(Antibody = gsub("^\\s*(.*$)", "\\1", Antibody),
                  Vendor = ifelse(Vendor == "BD Bioscience",
                                  "BD Biosciences", Vendor)) %>%
    dplyr::rename(Antigen = Antibody)

# Trzupek_2021 ----

trzupek_2021 <- adt_clones$Trzupek_2021 %>%
    dplyr::mutate(Study = "Trzupek_2021",
                  Cat_Number = as.character(Cat_Number))

# Ty 2023 ----
ty_2023 <- adt_clones$Ty_2023 %>%
    dplyr::rename(TotalSeq_Cat = TotalSeqCat,
                  Antigen = Antibody)

# Valenzi_2019 ----

# The only clone information provided is that the antibodies are TotalSeq A

# Vanuystel_2020 ----
vanuytsel <- adt_clones$Vanuytsel_2020 %>%
    dplyr::mutate(Study = "Vanuytsel_2020",
                  Cat_Number = as.character(Cat_Number)) %>%
    dplyr::select(-Gene)

# Wang_2020_PBMC ----

wang <- adt_clones$Wang_2020_PBMC
wang <- wang %>%
    dplyr::rename(Antigen = Antibody,
                  TotalSeq_Cat = Stain) %>%
    dplyr::mutate(Study = "Wang_2020_PBMC",
                  TotalSeq_Cat = gsub("TotalSeq", "", TotalSeq_Cat))

# Witkowski ----

witkowski <- adt_clones$Witkowski_2020
witkowski <- witkowski %>%
    dplyr::rename(Vendor = SOURCE,
                  temp = `REAGENT or RESOURCE`) %>%
    tidyr::separate(IDENTIFIER, into = c("Cat_Number", "Lot", "Clone"),
                    sep = "; Lot#|; Clone:") %>%
    dplyr::mutate(temp = gsub("CITE-Seq: ", "", temp),
                  Barcode_Sequence = gsub("^.*-([ACTG]+)$", "\\1", temp),
                  Antigen = stringr::str_replace(temp,
                                                 pattern = Barcode_Sequence,
                                                 replacement = ""),
                  Antigen = gsub("[ ]?-$", "", Antigen),
                  Cat_Number = gsub("Cat#", "", Cat_Number),
                  Study = "Witkowski_2020",
                  TotalSeq_Cat = "A") %>% # Determined by looking up BioLegend
    dplyr::filter(! grepl("Hashtag", temp)) %>%
    dplyr::select(-temp)
# Wu_2021 ----

wu_2021 <- adt_clones$Wu_2021
wu_2021 <- wu_2021 %>%
    dplyr::rename(Ensembl_ID = `Ensembl Gene Id`) %>%
    dplyr::mutate(Vendor = "BioLegend",
                  Study = "Wu_2021",
                  TotalSeq_Cat = gsub("^([ABC])\t?[0-9]+$", "\\1",
                                      Barcode_number),
                  Oligo_ID = gsub("^[ABC]\t?", "", Barcode_number),
                  Cat_Number = as.character(Cat_Number)) %>%
    # Biolegend_name is not necessary as Cat_Number is present
    dplyr::select(-Barcode_number, -Biolegend_name) %>%
    dplyr::select(-Panel_inclusion, -Barcode_V1, -Barcode_V2)


# 10x genomics data sets ----
tenx_id_to_fname <- tenx_id_to_fname %>%
    dplyr::mutate(filename = gsub("\\.csv", "", filename)) %>%
    dplyr::rename(Study = ID) %>% # TO DO: 10x19Nov2018-2 is inconsistent
    dplyr::select(Study, filename, Source, Notes)

tenx_clone_lns <- lapply(tenx_clones, nrow)
tenx_clones <- do.call(dplyr::bind_rows, tenx_clones) %>%
    dplyr::mutate(filename = rep(names(tenx_clones), tenx_clone_lns),
                  Vendor = "BioLegend",
                  TotalSeq_Cat = gsub(".*Total(Seq)?([ABCD])$", "\\2", name),
                  Clone = gsub(".*_([^_]*)_Total(Seq)?[ABCD]$", "\\1", name),
                  Clone = AbNames:::.noDups(Clone, name),
                  Clone = dplyr::na_if(Clone, "control")) %>%
    dplyr::rename(Antigen = id,
                  Barcode_Sequence = sequence,
                  TEMP = name) %>% # Renaming the data name as column
                  # "Data_Name" is removed prior to merging data
    dplyr::left_join(tenx_id_to_fname, by = "filename") %>%
    dplyr::select(-read, -pattern, -feature_type, -filename)

# Merge tables and make formatting more consistent ----

all_dat <- list(arunachalam, bai, buus, cadot, chung_2021, colpitts, fernandez,
                fidanza, frangieh, granja, hao, holmes_2020,
                kaufmann, kotliarov, lawlor, leader, lecoz, lee, li_2023, liu,
                mair, mimitou_2019, mimitou_2021, nathan, nettersheim, papalexi,
                pei_2020, poch_2021, pombo, pont, qi_2020, rincon,
                shangguan, stephenson, stoeckius_2017, stoeckius, stuart, su,
                triana_2021, trzupek_2020, trzupek_2021, ty_2023, vanuytsel,
                wang, witkowski, wu_2021, tenx_clones)

all_clones <- do.call(dplyr::bind_rows, all_dat) %>%
    dplyr::mutate(
        # remove unnecessary whitespace
        across(where(is.character), ~stringr::str_squish(.x)),

        Vendor = ifelse(Vendor == "Biolegend", "BioLegend", Vendor),

        Cat_Number = case_when(
            grepl("^(AHS)?[0-9-]+$", Cat_Number) ~ Cat_Number,
            # Common formatting error
            is.na(Cat_Number) ~ NA_character_,
            Cat_Number %in% c("n/a", "N/A") ~ NA_character_,
            Cat_Number == "N/A (Custom Order)" ~ "custom made",
            Cat_Number == "custom made" ~ "custom made",
            grepl("^\\*", Cat_Number) ~
                gsub("^\\*([0-9]+).*$","custom made (similar to \\1)",
                     Cat_Number),
            grepl("Biozol|R&D|Fluidigm", Vendor) ~ Cat_Number,
            TRUE ~ "ERROR"),

        RRID = ifelse(RRID %in% c("n/a", "N/A"), NA_character_, RRID),

        Reactivity = ifelse(Reactivity %in% c("", "n/a", "N/A"),
                            NA_character_, Reactivity),

        Reactivity = tolower(Reactivity),
        Reactivity = gsub(", ", "/", Reactivity),

        Gene_Name = na_if(Gene_Name, "/"),

        # Import errors:
        Clone = ifelse(Clone %in% c("50000000000", "5E+10"), "5E10", Clone),
        Clone = ifelse(Clone == "500000000", "5E8", Clone),
        Clone = ifelse(Clone == "9E+02", "9E2", Clone),
        Clone = ifelse(Clone %in% c("1E+03", "1000"), "10E2", Clone),

        Antigen = replaceGreekSyms(Antigen, 'sym2letter')) %>%

    # If Isotype column includes species, put the species into "Origin_Species"
    tidyr::separate(Isotype, into = c("ORI_SP", "Isotype"),
                    sep = " (?=Ig)", fill = "left") %>%
    dplyr::mutate(Origin_Species = dplyr::coalesce(Origin_Species, ORI_SP)) %>%
    dplyr::select(-ORI_SP)

# Check expected number of rows
sum(sapply(all_dat, nrow)) == nrow(all_clones)

# Manually identify controls ----

# Note: IgM, IgA, IgD and Ig light chain with reactivity against
# human are not isotype controls

all_clones <- all_clones %>%
    dplyr::mutate(Isotype_Control =
                    grepl("Mouse|Ra[tg]|[Hh]am|[Cc]ontrol|Ctrl|[Ii]so",
                          Antigen),

                  # Note this won't identify all of the isotype controls
                  Reactivity = case_when(
                    ! is.na(Reactivity) ~ Reactivity,
                    Isotype_Control & grepl("Mouse|^mIg", Antigen) ~ "mouse",
                    Isotype_Control & grepl("Rat|^rIg", Antigen) ~ "rat",
                    Isotype_Control &
                      grepl("Armenian Hamster|^[Aa]rm", Antigen) ~
                         "armenian hamster",
                    TRUE ~ Reactivity))

# Fix manually identified mistakes or inconsistencies ----

# Some of these are due to import errors, not in the original data

mistakes_clone_oligo <- tibble::tribble(
    ~Study, ~Antigen, ~Cat_Number, ~Clone, ~Oligo_ID, ~TotalSeq_Cat,
    "Arunachalam_2020", "FCER1a", "334641", "AER-37 (CRA-1)", "0352", "A",
    "Buus_2021", "CD103", "350233", "Ber-ACT8", "0145", "C",
    "Hao_2021", "CCR7", NA, "G043H7", NA, "C",
    "Mimitou_2019", "CCR7", NA, "G043H7", NA, NA,
    "Shangguan_2021","CD303","354241","201A","0370", "C",
    "Shangguan_2021","CD56","362559","5.1H11", NA, "C",
    "Witkowski_2020", "CD274", "329743", "29E.2A3", "0007", "A",
    "Witkowski_2020", "CD161", "339945", "HP-3G10", "0149", "A",
    "Witkowski_2020", "CD115", "347325", "9-4D2-1E4", "0398", "A",
    "Witkowski_2020", "CD134", "350033", "Ber-ACT35 (ACT35)", "0158", "A",
    "Witkowski_2020", "CD223", "369333", "11C3C65", "0152", "A",
    "Witkowski_2020", "CD152", "369619", "BNI3", "0151", "A",
    "Wu_2021", "CD45", "304064", "HI30", "0391", "A",
    "Wu_2021", "CD28", "302955", "CD28.2", "0386", "A",
    "Wu_2021", "CD127", "351352", "A019D5", "0390", "A",
    "Wu_2021", "CD10", "312231", "HI10a", "0062", "A",
    "Wu_2021", "OX40", "350033", "Ber-ACT35 (ACT35)", "0158", "A",
    "Wu_2021", "CD21", "354915", "Bu32", "0181", "A",
    "Wu_2021", "CD66b", "392905", "6/40c", "0166", "A",
    "Buus_2021", "CD134", "350035", "Ber-ACT35 (ACT35)", "0158", "C",
    "Qi_2020", "CD38", "940466", "HB-7", "AHS0189", NA,
    "Frangieh_2021", "CD274", NA, "29E.2A3", "0007", NA,
    "Leader_2021", "CD8", "301002", "RPA-T8", NA, NA,
    "PomboAntunes_2021", "CD178", "306413", "NOK-1", "0177", "A",
    "Liu_2021", "CD178 (Fas-L)", NA, "NOK-1", "0177", "C",
    # Actually an alternative name 2331 (FUN-1)
    "Trzupek_2021", "CD86", "940025", "2331", NA, NA
    )

mistakes_cat_num <- tibble::tribble(
    ~Study, ~Antigen, ~Cat_Number, ~Clone, ~Oligo_ID, ~TotalSeq_Cat,
    "Granja_2019", "CD4", "344649", "SK3", "0045", "A",
    "Granja_2019", "CD14", "301855","M5E2", "0081", "A",
    "Liu_2021", "anti-Tau Phospho (Thr181)", NA, "M7004D06", "0856", "C",
    "Liu_2021", "CD226 (DNAM-1)", "338337", "11A8", "0368", "C",
    "Liu_2021", "CD257 (BAFF, BLYS)", "366513", "1D6", "0863", "C",
    "Kotliarov_2020","CD314 (NKG2D)", "320835", "1D11", "0165", "A",
    "PomboAntunes_2021", "CD49d", "304337", "9F10", "0576", "A",
    "PomboAntunes_2021", "CD83", "305339", "HB15e", "0359", "A"
  )

mistakes_RRID <- tibble::tribble(
    ~Study, ~Antigen, ~Cat_Number, ~RRID,
    "Liu_2021", "CD226 (DNAM-1)", "338337", "AB_2800899",
    "Liu_2021", "CD257 (BAFF, BLYS)", "366513", "AB_2832723",
    "PomboAntunes_2021", "CD49d", "304337", "AB_2783166",
    "PomboAntunes_2021", "CD83", "305339", "AB_2800784"
)

all_clones <- all_clones %>%
    dplyr::rows_update(mistakes_clone_oligo,
                       by = c("Study", "Antigen", "Cat_Number",
                              "TotalSeq_Cat")) %>%
    dplyr::rows_update(mistakes_cat_num,
                       by = c("Study", "Antigen", "Clone",
                              "Oligo_ID", "TotalSeq_Cat")) %>%
    dplyr::rows_update(mistakes_RRID,
                       by = c("Study", "Antigen", "Cat_Number"))


# Fill in missing catalogue numbers using TotalSeq data set -----
# (note I am assuming that catalogue numbers are stable for
# Clone/Oligo/sequencing chemistry combination)

data(totalseq)
totalseq <- as_tibble(totalseq)

ts <- totalseq %>%
    dplyr::select(Clone, Oligo_ID, TotalSeq_Cat, Cat_Number) %>%
    unique()

all_clones <- all_clones %>%
    dplyr::rows_patch(ts, by = c("Clone", "Oligo_ID", "TotalSeq_Cat"),
                      unmatched = "ignore")

# Write table ----
# For the AbNames demonstration data set, we only want human antibodies

all_clones <- all_clones %>%
    dplyr::filter(grepl("[Hh]uman", Reactivity) |
                    is.na(Reactivity) |
                    Isotype_Control)

all_clones <- all_clones %>%
    dplyr::relocate(Study, Antigen, Clone, Cat_Number,
                    Oligo_ID, TotalSeq_Cat) %>%
    dplyr::arrange(Antigen, Cat_Number)


readr::write_delim(all_clones,
                   file = sprintf("%s/merged_adt_clones.tsv", clone_dir),
                   delim = ",")

# Clean up
rm(list = setdiff(ls(), c(existing, "existing", "all_clones")))

# Notes:
# MOPC273 Hao IgG2a is an error?
# CD8 - 146Nd Leader_2021 146Nd is a metal, cat number doesn't have the metal
