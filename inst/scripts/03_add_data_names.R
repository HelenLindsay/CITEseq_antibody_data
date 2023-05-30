# This script manually matches the Antibody names with the Antibody Derived Tag
# (ADT) names used in the data.  It requires a file of protein names, generated
# by the CITE-seq data fetching pipeline, see
# https://github.com/bernibra/CITE-wrangling

# To do: TEMP may be the data name, remove after merging
# Note: Wu_2021_a provided data on request


if(! grepl("vignettes", getwd())) { setwd("./vignettes") }

# Libraries ----

library("dplyr")
library("AbNames")
library("readxl")
library("readr")

# Get filenames by study ----

#protein_names <- list.files("inst/extdata/protein_names",
#                            full.names=TRUE, recursive=TRUE)
#protein_short_nms <- gsub(".*features_", "", protein_names)
#protein_names <- data.frame(long = protein_names, short = protein_short_nms)
#nm_to_study <- readr::read_delim("inst/extdata/study_to_name.csv")
#studies <- gsub(".*sce-objects\\/([^\\/]*)\\/.*", "\\1", protein_sce)
#protein_names <- lapply(protein_sce, function(x) rownames(read_rds(x)))
#names(protein_names) <- protein_sce
#protein_names <- split(protein_names, studies)


# Get a list of the unique names used for ADT expression in each study
# Here we only keep files with human data
all_protein_names <- read_rds("../inst/extdata/protein_names.rds")
studies <- names(all_protein_names)
protein_files <- unlist(all_protein_names, recursive = FALSE)
keep_files <- ! grepl("[Mm]ouse", names(protein_files))
filename_to_study <- rep(studies, lengths(all_protein_names))
protein_names <- split(protein_files[keep_files],
                    filename_to_study[keep_files])
ab_fnames <- split(names(protein_files)[keep_files],
                   filename_to_study[keep_files])
protein_names <- lapply(protein_names, function(x) unique(unlist(x)))

# Match study with data file names ----
acc_to_name <- read_delim("../inst/extdata/metadata/papers.csv") %>%
    dplyr::mutate(Accession =
                    dplyr::coalesce(Accession, gsub("_", "", Name))) %>%
    dplyr::filter(Accession %in% studies)
acc_to_study <- structure(acc_to_name$Accession, names = acc_to_name$Name)

# Load citeseq data ----

#citeseq <- readr::read_delim(
#    file = "../inst/extdata/ADT_clones/merged_adt_clones.tsv") %>%
#    # If rerunning, remove existing data names
#    dplyr::select(-any_of("Data_Name"))

# Check for missing studies ----
setdiff(studies, acc_to_name$Accession)
setdiff(citeseq$Study, acc_to_name$Name)

# Original colnames
original_colnames <- colnames(citeseq)

# Check which studies list the same Antigen more than once
# (doesn't include studies where antigen is listed as _1 _2 etc)
citeseq %>%
    dplyr::select(Study, Antigen, Clone, TotalSeq_Cat, Oligo_ID) %>%
    dplyr::group_by(Study) %>%
    unique() %>%
    dplyr::filter(n() != n_distinct(Antigen)) %>%
    dplyr::filter(AbNames:::.dups(Antigen)) %>%
    dplyr::arrange(Study, Antigen)  %>%
    data.frame()

duplicated_ab <- c("Granja_2019", "Hao_2021", "Leader_2021",
                   "Mimitou_2021", "Papalexi_2021", "PomboAntunes_2021")

# Arunachalam_2020 ----

# There appears to be an error in the data name -
# CD57--Ber-ACT8-TSA - Ber-ACT8 is the clone for CD103, not CD57

nm <- "Arunachalam_2020"

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

arunachalam_2020 <- citeseq %>% dplyr::filter(Study == nm)

data_names <- tibble(Data_Name = data_names) %>%
    tidyr::separate(Data_Name, into = c("Antigen", "Clone"),
                    sep = "--", remove = FALSE) %>%
    dplyr::mutate(Clone = gsub("-TSA", "", Clone)) %>%
    dplyr::rows_update(data.frame(Antigen = "CD274_PD-L1", Clone = "29E.2A3"),
                       by = "Antigen")

# Before joining, test for differences

# Antigen is shared but clone is different
x <- data_names %>%
    dplyr::semi_join(arunachalam_2020, by = "Antigen") %>%
    dplyr::anti_join(arunachalam_2020, by = "Clone")

arunachalam_2020 %>%
    dplyr::filter(Antigen %in% x$Antigen) %>%
    dplyr::select(Antigen, Clone) %>%
    dplyr::left_join(data_names, by = "Antigen")

# Clone is the same but antigen is different
y <- data_names %>%
    dplyr::semi_join(arunachalam_2020, by = "Clone") %>%
    dplyr::anti_join(arunachalam_2020, by = "Antigen")

arunachalam_2020 %>%
    dplyr::filter(Clone %in% y$Clone) %>%
    dplyr::select(Antigen, Clone) %>%
    dplyr::left_join(data_names, by = "Clone")

arunachalam_2020 <- arunachalam_2020 %>%
    AbNames:::left_join_any(data_names, cols = c("Antigen", "Clone"))

# Should be TRUE
! any(is.na(arunachalam_2020$Data_Name))
all(data_names$Data_Name %in% arunachalam_2020$Data_Name)

# Buus_2021 ----

nm <- "Buus_2021"
buus_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = gsub("Iso", "", Antigen))

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]


# Check that all are accounted for
setdiff(data_names, buus_2021$Data_Name)
setdiff(buus_2021$Data_Name, data_names)

# Cadot_2020 ----
# Isotype control is in clone table but not in data

nm <- "Cadot_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

cadot_2020 <- citeseq %>%
  dplyr::filter(Study == nm) %>%
  dplyr::mutate(Data_Name = gsub("\\(.*", "", Antigen),
                Data_Name = ifelse(Data_Name == "CD8a", "CD8", Data_Name),
                Data_Name = ifelse(Data_Name %in% data_names, Data_Name, NA))

all(data_names %in% cadot_2020$Data_Name)

# Chung_2021 ----

# Chang includes a single antibody profiled in HeLa cells and
# 4 antibodies profiled in mouse hippocampus.

nm <- "Chung_2021"

# Select the human antibodies:

#fnames <- ab_fnames[[acc_to_study[[nm]]]]
#data_names <- read_rds(fnames[! grepl("mouse", fnames)])
data_names <- protein_names[[acc_to_study[[nm]]]]

chung_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = Antigen)

setdiff(data_names, chung_2021$Antigen)
setdiff(chung_2021$Antigen, data_names)


# Fernandez_2019 -----

nm <- "Fernandez_2019"
print(nm)

#data_names <- unique(unlist(
#    lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

fernandez_2020 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = Antigen)

# Check
setdiff(data_names, fernandez_2020$Antigen)
setdiff(fernandez_2020$Antigen, data_names)

# Fidanza_2020 ----

nm <- "Fidanza_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("\\.TotalA", "", data_names))

fidanza_2020 <- citeseq %>% dplyr::filter(Study == nm)

all(data_names$Antigen %in% fidanza_2020$Antigen)
all(fidanza_2020$Antigen %in% data_names$Antigen)

fidanza_2020 <- dplyr::full_join(fidanza_2020, data_names)

# Frangieh_2021 ----

nm <- "Frangieh_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

frangieh_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = gsub("[- ]", "_", Antigen),
                  Data_Name = gsub("_[Ii]sotype.*", "", Data_Name),
                  Data_Name = gsub("(CD[0-9]+)([A-Z]$)", "\\1\\L\\2",
                                   Data_Name, perl = TRUE),
                  Data_Name = gsub("HLA_ABC", "HLA_A", Data_Name))

setdiff(frangieh_2021$Data_Name,  data_names)
setdiff(data_names, frangieh_2021$Data_Name)

# Granja_2019 ----

# Note: Granja et al used two different totalseq preparations,
# some of which have shared antibodies
# From supplementary tables T3:
# PBMC and BBMC used TotalSeq B preparation
# CD34 and MPAL used TotalSeq A preparation

# Note that we need the full patient table to know if MPAL are BBMC or PBMC

# Granja clone table is missing CD34 and CD38
# Isotype controls are not in the data table

nm <- "Granja_2019"
print(nm)

granja_nms <- ab_fnames[[acc_to_study[[nm]]]]
granja_samples <- gsub(".*scADT_(.*)\\.rds", "\\1", granja_nms)
granja_ts <- ifelse(grepl("CD34|MPAL", granja_samples), "A", "B")

# According to the data names, the same 16 panel is used in all
#data_names <- lapply(granja_nms, read_rds)
data_names <- all_protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = unlist(data_names),
                         TotalSeq_Cat = rep(granja_ts, lengths(data_names)),
                         Sample = rep(granja_samples, lengths(data_names))) %>%
    unique()

granja_2019 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = gsub(" .*", "", Antigen),
                  Data_Name = gsub("CD8a", "CD8A", Data_Name))

# Antigens in data missing from clone table (none)
dplyr::anti_join(data_names, granja_2019, by = c("Data_Name", "TotalSeq_Cat"))
# Antigens in clone table missing from data (isotype controls)
dplyr::anti_join(granja_2019, data_names, by = c("Data_Name", "TotalSeq_Cat"))

granja_2019 <- granja_2019 %>%
    dplyr::full_join(data_names, by = c("Data_Name", "TotalSeq_Cat"),
                     na_matches = "never") %>%
    tidyr::fill(Study)

# Hao_2021 ---- (ECCITE-SEQ MISSING)----

# CHECK: TWO DIFFERENT BARCODES FOR CD45RO AND CD45RA IN ECCITE-SEQ

# TotalSeq_C = ECCITE-Seq

# Hao dat193a is available through the Fred Hutch portal and through GEO
# The portal includes a meta data table relating antibodies to data names
# for the TotalSeq A (CITE-seq) study.

# Notch-1 and Notch-3 are listed in the clone table, Notch-1 and Notch-2 in data
# By matching catalogue number, Notch-3 (Antibody name) is correct

# CD193 and CD90 occur in the antibody table but not in the data.

nm <- "Hao_2021"
#data_names <- (unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- all_protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names)

# Some clones got converted to scientific notation somewhere along the line
hao_portal_patch <- tibble(Cat_Number =  c("310729", "352207",
                                           "331943", "328135"),
                           Clone = c("5E8", "10E2", "9E2", "5E10"))

portal_fname <- "../inst/extdata/protein_names_extra/Hao_portal_antibody_info.csv"
hao_portal <- read_delim(portal_fname) %>%
    dplyr::rename(Data_Name = "#protein",
                  Cat_Number = Catalog,
                  Antigen = Specificity) %>%
    dplyr::mutate(Cat_Number = as.character(Cat_Number),
                  Data_Name = gsub("_", "-", Data_Name)) %>%
    dplyr::select(Data_Name, Cat_Number, Clone) %>%
    dplyr::rows_update(hao_portal_patch, by = "Cat_Number")

# All of the data names from the portal are in the data
setdiff(data_names$Data_Name,  hao_portal$Data_Name)
setdiff(hao_portal$Data_Name, data_names$Data_Name)

hao_2021 <- citeseq %>%
    dplyr::filter(Study == nm)

anti_join(hao_2021, hao_portal, by = c("Cat_Number", "Clone")) %>%
    dplyr::filter(! Assay == "ECCITE-seq") %>%
    dplyr::select(Antigen, Cat_Number, Clone)

anti_join(hao_portal, hao_2021, by = c("Cat_Number", "Clone"))

hao_2021 <- hao_2021 %>%
    dplyr::left_join(hao_portal, by = c("Cat_Number", "Clone")) %>%
    dplyr::mutate(Data_Name = ifelse(Assay == "ECCITE-seq", NA, Data_Name))

# Holmes_2020 ----

nm <- "Holmes_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    dplyr::mutate(Antigen = gsub("(_Human)?\\..*", "", Data_Name),
                  Antigen = gsub(".*_([^_]+)", "\\1", Antigen))

holmes_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

# All of the data names from the portal are in the data
setdiff(data_names$Antigen,  holmes_2020$Antigen)
setdiff(holmes_2020$Antigen, data_names$Antigen)

holmes_2020 <- holmes_2020 %>%
  dplyr::left_join(data_names, by = "Antigen")

# Kaufmann_2021 ----

nm <- "Kaufmann_2021"
print(nm)

#data_names <- (unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    dplyr::mutate(Short = gsub("ADT-", "", Data_Name),
                Short = gsub("-", "", Short))

kaufmann_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Short = gsub("^(.*)\\s\\((.*)\\)$", "\\2-\\1", Antigen),
                Short = gsub("[-/]", "", Short),
                Short = gsub(".*integrin", "", Short),
                Short = AbNames::replaceGreekSyms(Short, "sym2letter"),
                Short = gsub("TCR ", "TCR", Short))

# Check that all antibodies are matched in Short
setdiff(data_names$Short, kaufmann_2021$Short)
setdiff(kaufmann_2021$Short, data_names$Short)

kaufmann_2021 <- kaufmann_2021 %>%
    dplyr::left_join(data_names, by = "Short") %>%
    dplyr::select(-Short)

# Kotliarov_2020 -----
# CD206 is present in data but missing from clone table

nm <- "Kotliarov_2020"
print(nm)

#data_names <- (unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("_PROT", "", data_names)) %>%
    dplyr::mutate(TEMP = stringr::str_squish(gsub("\\s", "", TEMP)),
                  TEMP = AbNames::replaceGreekSyms(TEMP, "word2letter"))

# Join in the data name, fill the necessary fields for the missing antigen
kotliarov_2020 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub("\\s\\(.*", "", Antigen),
                  TEMP = gsub("CX3CR1.*", "CX3CR1", TEMP),
                  TEMP = AbNames::replaceGreekSyms(TEMP, "word2letter"),
                  TEMP = gsub("\\s|,|(Ctrl)|\\/", "", TEMP))

setdiff(data_names$TEMP, kotliarov_2020$TEMP)
setdiff(kotliarov_2020$TEMP, data_names$TEMP)

kotliarov_2020 <- kotliarov_2020 %>%
    dplyr::full_join(data_names) %>%
    dplyr::mutate(Antigen = ifelse(is.na(Antigen), TEMP, Antigen)) %>%
    dplyr::select(-TEMP) %>%
    tidyr::fill(Study)


# Krebs 2020 ----
# Vendor info from paper supplement.
# Krebs is not in citeseq data set because clone info is not provided.

nm <- "Krebs_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

krebs_2020 <- data.frame(Antigen = gsub("ADT-", "", data_names),
                         Vendor = "BioLegend",
                         Data_Name = data_names,
                         Study = nm) %>%
    AbNames::addID()

# Lawlor_2021 ----
# CD185 is listed in the clone table but doesn't appear in the data

nm <- "Lawlor_2021"
print(nm)

#data_names <- unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("-[ACTG]+$", "", data_names))

lawlor_2021 <- citeseq %>%
    dplyr::filter(Study == nm)

setdiff(data_names$Antigen, lawlor_2021$Antigen)
setdiff(lawlor_2021$Antigen, data_names$Antigen)

lawlor_2021 <- lawlor_2021 %>%
    dplyr::left_join(data_names, by = "Antigen")

# Leader_2021 ---- TO DO change supp table location, check if TOTALSEQ filled previously----
# Leader uses a mixture of commercially conjugated (BioLegend) and
# in-house conjugated antibodies.  Different preparations of the same antibody
# are possible.

# From the supplementary data, we have panels containing antigen and barcode
# sequence for each patient.  We assume that when the barcode sequence does
# not match the expected BioLegend TotalSeq sequence, an antibody is an
# in-house conjugation.  We further assume that the same antibody panel was
# used when a patient sample was sequenced in several batches.

# I am unsure about matching these, e.g. EPCAM (= CD326) has a TotalSeq barcode
# in the patient table but CD326 catalogue number from the clone table is not
# a TotalSeq antibody

# Some of the samples have only RNA, so the protein name files are empty


nm <- "Leader_2021"

#leader_fnames <- ab_fnames[[acc_to_study[[nm]]]]
#data_names <- lapply(leader_fnames, read_rds)
data_names <- all_protein_names[[acc_to_study[[nm]]]]

# Get patient and batch ID from filenames
batch_id <- rep(gsub(".*batch_ID_([0-9]+).*", "\\1", names(data_names)),
                lengths(data_names))
patient_id <- rep(gsub(".*patient_([0-9]+).*", "\\1", names(data_names)),
                  lengths(data_names))
data_names <- data.frame(Data_Name = unlist(data_names)) %>%
    dplyr::mutate(batch = batch_id,
                  patient = patient_id)

# Check that antibodies are not duplicated within patient / batch
data_names %>%
    group_by(patient, Data_Name, batch) %>%
    arrange(patient, Data_Name) %>%
    filter(n() > 1)

# Load per-patient CITE-seq panels

cite_panels <- "../inst/extdata/protein_names_extra/Leader_2021_citeseq_panels.xlsx"

cite_sheets <- excel_sheets(cite_panels)
cite_panels <- lapply(cite_sheets[2:length(cite_sheets)], function(nm){
    readxl::read_excel(cite_panels, sheet = nm) %>%
        dplyr::mutate(patient = gsub("patient_", "", nm))
})

cite_panels <- dplyr::bind_rows(cite_panels) %>%
    dplyr::select(name, sequence, patient) %>%
    dplyr::arrange(name) %>%
    dplyr::rename(Antigen = name,
                  Barcode_Sequence = sequence)

# Add barcode sequence from patient CITE-seq panels into data_names table
data_names <- data_names %>%
    dplyr::left_join(cite_panels, by = c(Data_Name = "Antigen", "patient"))

# Check that one barcode is only ever associated with one data name
data_names %>%
    dplyr::group_by(Barcode_Sequence) %>%
    dplyr::filter(n_distinct(Data_Name) > 1)

# Patient and batch are no longer necessary
data_names <- data_names %>%
    dplyr::select(-patient, -batch) %>%
    unique()


# Add the barcode sequence into the clone table using the TotalSeq vendor
# information
data(totalseq, package = "AbNames")
totalseq <- totalseq %>%
    dplyr::select(Cat_Number, Barcode_Sequence)

leader_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    # Barcode sequence is empty, add it from totalseq
    dplyr::select(-Barcode_Sequence) %>%
    dplyr::left_join(totalseq, by = "Cat_Number")

# Check that all TotalSeq antibodies have a barcode sequence
leader_2021 %>% filter(! is.na(TotalSeq_Cat) & is.na(Barcode_Sequence))

# Match clone to data name using barcode sequence
leader_2021 <- leader_2021 %>%
    dplyr::left_join(data_names, by = "Barcode_Sequence")

# Assume antigen matches if Barcode_Sequence is NA and Antigen matches Data_Name
# Select unmatched via barcode sequence in data table
data_names_patch <- data_names %>%
    dplyr::anti_join(leader_2021, by = "Barcode_Sequence") %>%
    # Assume that the short barcodes are the in-house conjugations
    dplyr::filter(nchar(Barcode_Sequence) == 6)

leader_2021 <- leader_2021 %>%
    dplyr::rows_patch(data_names_patch %>% dplyr::rename(Antigen = Data_Name),
                      by = c("Antigen"), unmatched = "ignore") %>%
    dplyr::rows_patch(data_names_patch,
                      by = "Barcode_Sequence", unmatched = "ignore")

# Unmatched via barcode sequence in clone table
data_no_match <- data_names %>%
  dplyr::anti_join(leader_2021, by = "Barcode_Sequence")

# Join in the remaining data names
leader_2021 <- leader_2021 %>%
    dplyr::full_join(data_no_match) %>%
    tidyr::fill(Study) %>%
    dplyr::mutate(Antigen = dplyr::coalesce(Antigen, Data_Name))


# LeCoz_2021 ----
# According to the paper, there are 8 isotype controls
# These are not included in the data
# The TotalSeq cocktail from the BioLegend website has 7 controls plus CD44

nm <- "LeCoz_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("anti-[a-z-]+_", "", data_names)) %>%
    dplyr::mutate(TEMP = gsub("_", " ", TEMP),
                  TEMP = ifelse(grepl("^CD[0-9]+[Pabe1]* ", TEMP),
                                gsub(" .*", "", TEMP),
                                TEMP),
                  TEMP = gsub("[ -]", "",
                              replaceGreekSyms(TEMP, "word2letter")))

lecoz_2021 <- citeseq %>%
  dplyr::filter(Study == nm) %>%
  dplyr::mutate(TEMP = ifelse(grepl("^CD[0-9]+[Pabe1]* ", Antigen),
                                    gsub(" .*", "", Antigen),
                                    Antigen),
                TEMP = gsub("[-,\\/ ]", "", TEMP))

setdiff(data_names$TEMP, lecoz_2021$TEMP)
setdiff(lecoz_2021$TEMP, data_names$TEMP)

lecoz_2021 <- lecoz_2021 %>%
    dplyr::full_join(data_names, by = "TEMP") %>%
    dplyr::select(-TEMP)

# Lee_2021 ----

nm <- "Lee_2021"
#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("-[-]?.*", "", data_names))

lee_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub("[ -].*", "", Antigen))

setdiff(data_names$TEMP, lee_2021$TEMP)
setdiff(lee_2021$TEMP, data_names$TEMP)

lee_2021 <- lee_2021 %>%
    dplyr::full_join(data_names, by = "TEMP") %>%
    dplyr::select(-TEMP)

# Liu_2021 ----

# I have matched (data) "TCR" with (clone table) "TCR Vbeta13.1" by elimination,
# but I'm not completely sure about this.

nm <- "Liu_2021"
#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

# Check that it's safe to strip off the .1 ending from data names
# any(AbNames:::.dups(gsub("\\.1", "", data_names)))

clone2data = tibble::tribble(~TEMP, ~Antigen,
                             "c-Met", "anti-c-Met",
                             "Phospho-Tau", "anti-Tau",
                             "IgG1Kiso", "Mouse IgG1 kappa isotype Ctrl",
                             "IgG2aKiso", "Mouse IgG2a kappa isotype Ctrl",
                             "IgG2bKiso", "Mouse IgG2b kappa isotype Ctrl",
                             "IgGFc", "IgG Fc",
                             "ratIgG2bKiso", "Rat IgG2b kappa Isotype Ctrl",
                             "IgKappaLight", "Ig light chain kappa",
                             "IgLambdaLight", "Ig light chain lambda",
                             "CD79b", "CD79b (IgBeta)",
                             "TCRab", "TCR alphabeta",
                             "TCRgd", "TCR gammadelta",
                             "TCRVa24-Jz18", "TCR Valpha24-Jalpha18 (iNKT cell)",
                             "TCRVa7.2", "TCR Va7.2",
                             "TCRVd2", "TCR Vdelta2",
                             "TCR", "TCR Vbeta13.1",
                             "TCRVg9", "TCR Vgamma9")

data_names <- data.frame(Data_Name = data_names[!grepl("^HTO", data_names)]) %>%
    dplyr::mutate(TEMP = gsub("\\.1", "", Data_Name)) %>%
    dplyr::full_join(clone2data, by = "TEMP") %>%
    dplyr::mutate(TEMP = dplyr::coalesce(Antigen, TEMP)) %>%
    dplyr::select(-Antigen)

liu_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = ifelse(grepl("Ig|TCR", Antigen),
                                Antigen,
                                gsub(" .*", "", Antigen)),
                  TEMP = gsub("[,\\/]", "", TEMP))

setdiff(data_names$TEMP, liu_2021$TEMP)
setdiff(liu_2021$TEMP, data_names$TEMP)

liu_2021 <- liu_2021 %>%
    dplyr::full_join(data_names, by = "TEMP") %>%
    dplyr::select(-TEMP)


# Mair_2020 ----
# CD223 in data names is Lag3 in clone table
# CD184 is in data not clone table

nm <- "Mair_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    dplyr::filter(grepl("pAbO$", Data_Name)) %>%
    tidyr::separate(Data_Name, into = c("Antigen", "Clone", "RRID", "Tag"),
                    remove = FALSE, sep = "\\|") %>%
    dplyr::mutate(Antigen = ifelse(Antigen == "CD223", "Lag3", Antigen))

mair_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

setdiff(data_names$Antigen, mair_2020$Antigen)
setdiff(mair_2020$Antigen, data_names$Antigen)

mair_2020 <- mair_2020 %>%
    dplyr::full_join(data_names %>% dplyr::select(Data_Name, Antigen),
                     by = "Antigen") %>%
    tidyr::fill(Study)

# Mimitou_2019 ----
#
# These antibodies appear in the K562 CRISPR screen.
# The barcodes match the antibodies used in the species mix experiment,
# where they are reported as hCD29 and hCD46.
# "CD29.AATAGCGGAGCC", "CD46.GCCAATTGCACT"

# This ADT is present in data but missing from clone table:
# "hCD13" (species mix experiment)

nm <- "Mimitou_2019"
print(nm)

data_names_patch <- tibble::tribble(~Data_Name, ~TEMP,
                                    "CD29.AATAGCGGAGCC", "hCD29",
                                    "CD46.GCCAATTGCACT", "hCD46")

mimitou_patch <- tibble::tribble(~Antigen, ~TEMP,
                                 "CD3E", "CD3",
                                 "CD8A", "CD8",
                                 "PD-L1", "B7.H1..PD.L1.")

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = data_names) %>%
    dplyr::rows_update(data_names_patch, by = "Data_Name")

mimitou_2019 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub("(.*) \\(mouse\\)", "m\\1", Antigen),
                  TEMP = gsub("(.*) \\(human\\)", "h\\1", TEMP),
                  TEMP = gsub("[^[:alnum:]]", "\\.", TEMP)) %>%
    dplyr::rows_update(mimitou_patch, by = "Antigen")

setdiff(data_names$TEMP, mimitou_2019$TEMP)
setdiff(mimitou_2019$TEMP, data_names$TEMP)

# Note: 2 rows gained because same antibody has different data names,
# 1 row gained for antibody not present in clone table
mimitou_2019 <- mimitou_2019 %>%
    dplyr::full_join(data_names, by = "TEMP") %>%
    tidyr::fill(Study) %>%
    dplyr::select(-TEMP) %>%
    dplyr::mutate(Antigen = ifelse(is.na(Antigen), "CD13", Antigen))


# Mimitou_2021 NOT ABLE TO MATCH ANTIGEN WITH DIFFERENT CLONES ----

# I emailed to request the barcode sequences for distinguishing clones
# but have not received an answer.

mimitou_patch <- tibble::tribble(~Data_Name, ~TEMP,
                                 "Cadherin11", "Cadherin 11",
                                 # Assuming these are the same in
                                 # different subexperiments
                                 "CD56.NCAM.Recombinant", "CD56 Recombinant",
                                 "CD56(NCAM)Recombinant", "CD56 Recombinant",

                                 "LOX.1", "LOX-1",
                                 "CD57_Recombinant", "CD57 Recombinant",
                                 "Folate_Receptor", "Folate Receptor b",
                                 "integrin_7", "integrin b7",
                                 "Notch1", "Notch 1",
                                 "Notch3", "Notch 3",
                                 "Rat_IgG1k_isotypeCtrl",
                                     "Rat IgG1, k isotype Ctrl",
                                 "Rat_IgG1l_IsotypeCtrl",
                                     "Rat IgG1, l Isotype Ctrl",
                                 "Rat_IgG2b_IsotypeCtrl",
                                     "Rat IgG2b, k Isotype Ctrl",
                                 "Rat_IgG2c_IsotypeCtrl",
                                     "Rat IgG2c, k Isotype Ctrl",
                                 "TCR_V_24-J_18", "TCR Va24-Ja18",
                                 "TCR_V_24.J_18", "TCR Va24-Ja18",
                                 "TCR_V_9", "TCR Vg9",
                                 "TCR_V_7.2", "TCR Va7.2",
                                 "TCRab", "TCR a/b",
                                 "TCRgd", "TCR g/d",
                                 "TCR_V_2", "TCR Vd2",
                                 "TIM.4", "TIM-4",
                                 "VEGFR3", "VEGFR-3")


nm <- "Mimitou_2021"
#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("\\(.*\\)", "", data_names)) %>%
    dplyr::rows_update(mimitou_patch, by = "Data_Name")

mimitou_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub(" ?\\(.*\\)", "", Antigen)) %>%
    # If there is no match because of multiple clones,
    # add entries without clone info
    dplyr::full_join(data_names, by = "TEMP") %>%
    tidyr::fill(Study) %>%
    dplyr::mutate(TEMP = gsub("[\\.-][12]$", "", TEMP),
                  TEMP = gsub("\\.", " ", TEMP),
                  Antigen = dplyr::coalesce(Antigen, TEMP)) %>%
    dplyr::select(-TEMP)


# Nathan_2021 ----

nm <- "Nathan_2021"
print(nm)

nathan_patch <- data.frame(Data_Name = c("MouseIgG_protein", "TCRab_protein"),
                           TEMP = c("Mouse", "TCR"))

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("[_\\.].*", "", data_names)) %>%
    dplyr::rows_update(nathan_patch, by = "Data_Name")

nathan_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub("[ \\/-].*", "", Antigen))

setdiff(data_names$TEMP, nathan_2021$TEMP)
setdiff(nathan_2021$TEMP, data_names$TEMP)

nathan_2021 <- nathan_2021 %>%
    dplyr::full_join(data_names, by = "TEMP") %>%
    dplyr::select(-TEMP)

# Papalexi_2021 ----

# Some of the antibodies in the data are not listed in the clone table

nm <- "Papalexi_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names)

papalexi_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = gsub("-", "", Antigen))

setdiff(data_names$Data_Name, papalexi_2021$Data_Name)
setdiff(papalexi_2021$Data_Name, data_names$Data_Name)

papalexi_2021 <- papalexi_2021 %>%
    dplyr::full_join(data_names, by = "Data_Name") %>%
    tidyr::fill(Study)

# Pei_2020 ----

# Some of the CD antigens have a 1 added at the end of the data name
# (all of the non-identical entries match this pattern,
# none of the antibodies are repeated)

nm <- "Pei_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("_ADT", "", data_names))

pei_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

suffix_1 <- setdiff(data_names$Antigen, pei_2020$Antigen)
data_names <- data_names %>%
    dplyr::mutate(Antigen = ifelse(Antigen %in% suffix_1,
                                   gsub("1$", "", Antigen), Antigen),
                  Antigen = gsub("\\.", "-", Antigen))

setdiff(data_names$Antigen, pei_2020$Antigen)
setdiff(pei_2020$Antigen, data_names$Antigen)

pei_2020 <- dplyr::left_join(pei_2020, data_names, by = "Antigen")

# Poch_2021 ----
# The data names are usually the HGNC symbol with a few exceptions
# There are some antibodies in the data missing from the clone table

nm <- "Poch_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Study = "Poch_2021",
                         Antigen = data_names)

poch_patch <- data.frame(Antigen = c("HLA-DR", "TCRgd", "CD8"),
                         Data_Name = c("HLA", "TCRgd", "CD8"))

poch_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = HGNC_SYMBOL) %>%
    dplyr::rows_update(poch_patch, by = "Antigen")

poch_aj <- dplyr::anti_join(data_names, poch_2021, by = "Data_Name")
poch_2021 <- poch_2021 %>%
    dplyr::bind_rows(poch_aj)

# PomboAntunes_2021 ----

# Study includes a mouse and and a human panel, here only matching the human
# Includes the same antibody with different oligo ID
# Isotype controls are not in the human data.
# Integrin beta 7 is listed as Integrin-b7 (mouse) IntegrinBeta7 (human) -
# (based on capitalisation of gene name Itgb7 / ITGB7)
# Names in Reactivity indicate panel e.g. mouse/human versus human/mouse
# Some isotype controls appear in the clone table but not the data,
# possibly used in the mouse panel

pombo_patch <- data.frame(Data_Name = c("B7.H4",
                                        "CD11a.CD18",
                                        "CD235a.CD235b",
                                        "CD307c",
                                        "DopamineReceptorD4",
                                        "Ig.lightChainK",
                                        # Below was originally Ig.lightChainλ
                                        "Ig.lightChain..",
                                        "LOX.1",
                                        "mastCellTryptase",
                                        "Tim.4",
                                        "IntegrinB7",
                                        "TCR.Va24.Ja18",
                                        "TCR.Va7.2",
                                        "TCR.Vd2",
                                        "TCR.Vg9"),
                          TEMP = c("B7H4",
                                   "CD11a_CD18 (LFA_1)",
                                   "CD235ab",
                                   "CD307c_FcRL3",
                                   "DopamineD4receptor",
                                   "Ig_light_chain_k",
                                   "Ig_light_chain_l",
                                   "LOX1",
                                   "Mast Cell Tryptase",
                                   "Tim4",
                                   "IntegrinBeta7",
                                   "TCRa24",
                                   "TCRa7",
                                   "TCRd2",
                                   "TCRg9"))

nm <- "PomboAntunes_2021"

#pombo_fnames <- ab_fnames[[acc_to_study[[nm]]]]
#"GSE163120_GSE163120_Human-GSM4972212_Citeseq_Human"
#pombo_human <- pombo_fnames[grepl("Human", pombo_fnames)]
pombo_human <- all_protein_names[[acc_to_study[[nm]]]]
pombo_human <- pombo_human[grepl("Human", names(pombo_human))]

# Oligo ID is present if there are two different clones used
data_names <- data.frame(Data_Name = unname(unlist(pombo_human))) %>%
    dplyr::mutate(Oligo_ID = gsub(".*\\.A([0-9]{4})$", "\\1", Data_Name),
                  Oligo_ID = AbNames:::.noDups(Oligo_ID, Data_Name),
                  TEMP = gsub("\\.A[0-9]{4}$", "", Data_Name),
                  TEMP = gsub("[\\.+-]", "_", TEMP)) %>%
    dplyr::rows_update(pombo_patch, by = "Data_Name")

pomboAntunes_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub("[\\.+-]", "_", Antigen)) %>%
    dplyr::mutate(TEMP = case_when(grepl("CD207", Antigen) &
                Reactivity == "human" ~ "CD207_h",
                                   grepl("CD207", Antigen) &
                Reactivity %in% c("mouse/human", "human/mouse") ~ "CD207_mh",
                                   TRUE ~ TEMP))

# Everything in the data is now in TEMP
x <- pomboAntunes_2021 %>%
    dplyr::filter(grepl("^human|isotype", Reactivity) &
                    (Gene_Name == toupper(Gene_Name) | is.na(Gene_Name)) ) %>%
    dplyr::pull(TEMP)

setdiff(x, data_names$TEMP)
setdiff(data_names$TEMP, x)

pomboAntunes_2021 <- pomboAntunes_2021 %>%
    # First add data names that have Oligo_IDs
    dplyr::left_join(data_names, by = c("TEMP", "Oligo_ID"),
                     na_matches = "never") %>%
    # Then add data names without Oligo_IDs
    dplyr::rows_patch(data_names %>% dplyr::filter(is.na(Oligo_ID)),
                     by = "TEMP", unmatched = "ignore") %>% # FIX CD207, REMOVE
    dplyr::select(-TEMP) %>%
    tidyr::fill(Study)

# Consider adding the mouse samples in future.
# pombo_mouse <- pombo_fnames[grepl("Mouse", pombo_fnames)]

# Pont_2020 -----

nm <- "Pont_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

pont_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("_Ab", "", data_names)) %>%
    dplyr::mutate(Antigen = ifelse(grepl("Isoytpe", Antigen),
                                   gsub("Isoytpe m(.*)",
                                        "Mouse \\1, kappa isotype Ctrl",
                                        Antigen), Antigen))
pont_2020 <- citeseq %>%
  dplyr::filter(Study == "Pont_2020")

setdiff(data_names$Antigen, pont_2020$Antigen)
setdiff(pont_2020$Antigen, data_names$Antigen)

pont_2020 <- dplyr::full_join(pont_2020, data_names)

# Qi_2020 ----

# Check it's safe to remove .1 suffix
#any(AbNames:::.dups(gsub("\\.1", "", data_names)))

nm <- "Qi_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("CD11b/Mac-1", "CD11b",
                                        gsub("\\.1", "", data_names)))

qi_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

setdiff(data_names$Antigen, qi_2020$Antigen)
setdiff(qi_2020$Antigen, data_names$Antigen)

qi_2020 <- qi_2020 %>%
    dplyr::full_join(data_names, by = "Antigen")

# RinconArevalo_2021 ----

nm <- "RinconArevalo_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

rinconArevalo_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = ifelse(Antigen == "CD4", "CD4.1", Antigen))

setdiff(data_names, rinconArevalo_2021$Data_Name)
setdiff(rinconArevalo_2021$Data_Name, data_names)

# Shangguan_2021 ----
# Isotype controls are not included in the data

nm <- "Shangguan_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

shangguan_2021 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(Data_Name = gsub(" ", "-", Antigen),
                  Data_Name = gsub("\\/|-\\(Fc\\)", "", Data_Name),
                  Data_Name = paste0(Data_Name, "-PROT"))

setdiff(data_names, shangguan_2021$Data_Name)
setdiff(shangguan_2021$Data_Name, data_names)

shangguan_2021 <- shangguan_2021 %>%
    dplyr::mutate(Data_Name = ifelse(Data_Name %in% data_names, Data_Name, NA))

# Stephenson_2021 ----
# Note that I have guessed TCR = TCR a/b and TCR_Vg2 = TCR Vd2

stephenson_patch <- data.frame(Data_Name = c("HLA-ABC",
                                             "Mouse IgG1_K_Iso",
                                             "Mouse_IgG2a_K_Iso",
                                             "Mouse_IgG2b_K_Iso",
                                             "Rat_IgG2b_K_Iso",
                                             "podoplanin",
                                             "TCRg_d",
                                             "FcERIa",
                                             "IgG_Fc",
                                             "TCR_Vg9",
                                             "HLA-A_2",
                                             "CEACAM1/5/6",
                                             "TCR_VB_13_1",
                                             "TCR_Va24-Ja18",
                                             "c-Met",
                                             "Iglambda",
                                             "Igkappa",
                                             "TCR_Va7.2",
                                             "KIR2DL5A",
                                             "FCGR2A",
                                             "phosphoTau",
                                             "GITR",
                                             "RANKL",
                                             "BAFF",
                                             "LIGHT",
                                             "BAFFR",
                                             "TACI",
                                             "MMR",
                                             "langerin",
                                             "OX40L",
                                             "PD1",
                                             "TCR_Vg2",
                                             "TCR"),
                        TEMP = c("HLA-A,B,C",
                                 "Mouse IgG1, k isotype Ctrl",
                                 "Mouse IgG2a, k isotype Ctrl",
                                 "Mouse IgG2b, k isotype Ctrl",
                                 "Rat IgG2b, k Isotype Ctrl",
                                 "Podoplanin",
                                 "TCR g/d",
                                 "FceRIa",
                                 "IgG Fc",
                                 "TCR Vg9",
                                 "HLA-A2",
                                 "CD66a/c/e",
                                 "TCR Vb13.1",
                                 "TCR Va24-Ja18",
                                 "anti-c-Met",
                                 "Ig light chain l",
                                 "Ig light chain k",
                                 "TCR Va7.2",
                                 "CD158f",
                                 "CD32",
                                 "Tau Phospho",
                                 "CD357",
                                 "CD254",
                                 "CD257",
                                 "CD258",
                                 "CD268",
                                 "CD267",
                                 "CD206",
                                 "CD207",
                                 "CD252",
                                 "CD279",
                                 "TCR Vd2",
                                 "TCR a/b")) %>%
    dplyr::mutate(Data_Name = sprintf("AB_%s", Data_Name))


nm <- "Stephenson_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("AB_", "", data_names)) %>%
    dplyr::rows_update(stephenson_patch)

stephenson_2021 <- citeseq %>%
   dplyr::filter(Study == nm) %>%
   dplyr::mutate(TEMP = gsub(" \\(.*", "", Antigen),
                 TEMP = ifelse(HGNC_SYMBOL %in% data_names$TEMP,
                               HGNC_SYMBOL, TEMP))

setdiff(data_names$TEMP, stephenson_2021$TEMP)
setdiff(stephenson_2021$TEMP, data_names$TEMP)

stephenson_2021 <- stephenson_2021 %>%
    dplyr::left_join(data_names, by = "TEMP") %>%
    dplyr::select(-TEMP)

# Stoeckius_2017 ----

nm <- "Stoeckius_2017"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

# CD29 is only in the mouse experiment and isn't in this data set
stoeckius_2017_patch <- data.frame(Antigen = c("CD3e", "CD8a","CD29"),
                                   Data_Name = c("CD3", "CD8", NA))

stoeckius_2017 <- citeseq %>%
    dplyr::filter(Study == nm)

setdiff(data_names, stoeckius_2017$Antigen)
setdiff(stoeckius_2017$Antigen, data_names)

stoeckius_2017 <- stoeckius_2017 %>%
    dplyr::mutate(Data_Name = Antigen) %>%
    dplyr::rows_update(stoeckius_2017_patch, by = "Antigen")

# Stoeckius_2018 ----
# Isotype controls are not in clone table

nm <- "Stoeckius_2018"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <-  data.frame(Data_Name = data_names,
                          Antigen = gsub("\\.[ACTG]+$", "", data_names)) %>%
    dplyr::mutate(Antigen = gsub("\\.", "-", Antigen))

stoeckius_2018 <- citeseq %>%
  dplyr::filter(Study == nm)

setdiff(stoeckius_2018$Antigen, data_names$Antigen)
setdiff(data_names$Antigen, stoeckius_2018$Antigen)

stoeckius_2018 <-
    dplyr::full_join(stoeckius_2018, data_names) %>%
    tidyr::fill(Study)

# Stuart_2019 ----

nm <- "Stuart_2019"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    dplyr::mutate(Antigen = ifelse(grepl("HLA", Data_Name),
                                   gsub("\\.", "-", Data_Name),
                                   gsub("\\..*", "", Data_Name)),
                  Antigen = gsub("CD278", "ICOS", Antigen))

stuart_2019 <- citeseq %>%
  dplyr::filter(Study == nm)

# Check that all are accounted for
setdiff(data_names$Antigen, stuart_2019$Antigen)
setdiff(stuart_2019$Antigen, data_names$Antigen)

# Merge in data names
stuart_2019 <- stuart_2019 %>%
    dplyr::left_join(data_names, by = "Antigen")

# Su_2020 ----

nm <- "Su_2020"
print(nm)

data2clone <- tibble::tribble(~Data_Name, ~REPLACE,
                              "CD158B_KIR2DL2-L3_NKAT3",
                              "CD158BKIR2DL2L3NKAT2",
                              "CD79_IgB", "CD79BIGB",
                              "ITGB7-1", "INTEGRINB7",
                              "LOX-1", "LOX1",
                              "Mouse_IgG1_Isotype-Control",
                              "MOUSEIGG1KISOTYPECTRL",
                              "Mouse_IgG2a_Isotype-Control",
                              "MOUSEIGG2AKISOTYPECTRL",
                              "Mouse_IgG2b_Isotype-Control",
                              "MOUSEIGG2BKISOTYPECTRL",
                              "Rat_IgG2b_Isotype-Control",
                              "RATIGG2BKISOTYPECTRL",
                              "TSLPR_CRL2", "TSLPRTSLPR")

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = toupper(data_names)) %>%
    # Strip off trailing -1 (there are no duplicated antibodies)
    # strip off punctuation and convert everything to upper case
    dplyr::mutate(TEMP = gsub("-1$", "", TEMP),
                  TEMP = gsub("[^[:alnum:]]", "", TEMP)) %>%
    dplyr::full_join(data2clone) %>%
    dplyr::mutate(TEMP = dplyr::coalesce(REPLACE, TEMP)) %>%
    dplyr::select(-REPLACE)


su_2020 <- citeseq %>%
  dplyr::filter(Study == nm) %>%
  dplyr::mutate(TEMP = AbNames::replaceGreekSyms(Antigen, "sym2letter"),
                TEMP = gsub("-1\\)$", "", TEMP),
                TEMP = gsub("[^[:alnum:]]", "", toupper(TEMP)))

setdiff(su_2020$TEMP, data_names$TEMP)
setdiff(data_names$TEMP, su_2020$TEMP)

su_2020 <- su_2020 %>%
  dplyr::full_join(data_names) %>%
  dplyr::select(-TEMP)

# Trzupek_2020 ----
# Oligo_ID was not in the clone table but can be filled in
# using the data names

nm <- "Trzupek_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    tidyr::separate(Data_Name, sep = "\\|", remove = FALSE,
                    into = c("TEMP", "HGNC_SYMBOL",
                             "Oligo_ID", "DELME")) %>%
    dplyr::mutate(HGNC_SYMBOL = gsub("_", ", ", HGNC_SYMBOL),
                  # Don't want to match using the HGNC_SYMBOL for CD45 isoforms
                  HGNC_SYMBOL = ifelse(grepl("CD45", TEMP), NA, HGNC_SYMBOL))

trzupek_2020 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = gsub(" .*", "", Antigen)) %>%
    dplyr::select(-Oligo_ID) %>% # Because it's NA in the clone table
    # Join by either HGNC_SYMBOL or TEMP name
    AbNames:::left_join_any(data_names, cols = c("HGNC_SYMBOL", "TEMP")) %>%
    dplyr::select(-DELME, -TEMP)

# Check nothing missing
setdiff(data_names$Data_Name, trzupek_2020$Data_Name)

# Trzupek_2021 ----

# Oligo_ID was not always in the clone table, we can fill in in using the
# data names
nm <- "Trzupek_2021"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[["Trzupek2021"]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    tidyr::separate(Data_Name, sep = "\\|", remove = FALSE,
                    into = c("TEMP", "HGNC_ID", "Oligo_ID", "delme")) %>%
    tidyr::separate(TEMP, into = c("TEMP", "Clone"),
                    sep =":", fill = "right") %>%
    dplyr::mutate(TEMP = ifelse(TEMP == "Tim3", "TIM3", TEMP))

trzupek_2021 <- citeseq %>%
    dplyr::filter(Study == "Trzupek_2021") %>%
    dplyr::mutate(TEMP = gsub(" .*", "", Antigen)) %>%
    dplyr::mutate(TEMP2 = gsub(".*\\((.*)\\)$", "\\1", Antigen),
                  TEMP = ifelse(TEMP %in% data_names$TEMP, TEMP, TEMP2)) %>%
    dplyr::select(-TEMP2)

setdiff(data_names$TEMP, trzupek_2021$TEMP)
setdiff(trzupek_2021$TEMP, data_names$TEMP)

trzupek_2021 <- trzupek_2021 %>%
    left_join(data_names %>% select(Data_Name, TEMP),
              by = c("TEMP")) %>%
    dplyr::rows_patch(data_names %>% select(Data_Name, Oligo_ID)) %>%
    dplyr::select(-TEMP) %>%
    dplyr::mutate(Study = "Trzupek_2021")


# Valenzi_2019 ----
# The only information given in the paper was that antibodies are TotalSeq A

nm <- "Valenzi_2019"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

valenzi_2019 <- data.frame(Antigen = gsub("-CITE", "", data_names),
                           Data_Name = data_names,
                           TotalSeq_Cat = "A",
                           Vendor = "BioLegend",
                           Study = nm)

setdiff(data_names, valenzi_2019$Data_Name)
setdiff(valenzi_2019$Data_Name, data_names)

# Vanuystel_2020 ----

vanuytsel_patch <- data.frame(Data_Name = c("VECAD-TCCACTCATTCTG",
                                            "TEK-protein-CGATCCCTTACCT",
                                            "FLT4-TGATCCGAAGTCG",
                                            "HLA-DR-DP-DQ-AATAGCGAGCAAG",
                                            "Pan-HLA-I-TATGCGAGGCTTA",
                                            "ENG-ATCGTCGAGAGCT"),
                              Antigen = c("CD144", "CD202b", "VEGFR3",
                                          "HLA-DR,DP,DQ", "pan-HLA class I",
                                          "CD105"))

nm <- "Vanuytsel_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("-protein", "", data_names)) %>%
    dplyr::filter(! Data_Name == "unmapped") %>%
    tidyr::separate(TEMP, into = c("Antigen", "Barcode_Sequence"),
                    sep = "-(?=[ACTG]+)") %>%
    dplyr::rows_update(vanuytsel_patch, by = "Data_Name")

vanuytsel_2020 <- citeseq %>%
    dplyr::filter(Study == nm)

setdiff(vanuytsel_2020$Antigen, data_names$Antigen)
setdiff(data_names$Antigen, vanuytsel_2020$Antigen)

vanuytsel_2020 <- vanuytsel_2020 %>%
    dplyr::full_join(data_names %>% select(-Barcode_Sequence),
                     by = "Antigen") %>%
    dplyr::rows_patch(data_names, by = c("Antigen", "Data_Name"))

# Wang_2020_PBMC ----

nm <- "Wang_2020_PBMC"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         Antigen = gsub("-[ACTG]+$", "", data_names))

wang_2020 <- citeseq %>%
  dplyr::filter(Study == nm)

# Check all accounted for
setdiff(data_names$Antigen, wang_2020$Antigen)
setdiff(wang_2020$Antigen, data_names$Antigen)

wang_2020 <- wang_2020 %>%
    dplyr::left_join(data_names, by = "Antigen")

# Witkowski_2020 ----
# Several present in data but missing from clone table.
# Data has ADT-PD1L1-CD274, clone has CD274 = PD-L1

nm <- "Witkowski_2020"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[[acc_to_study[[nm]]]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names,
                         TEMP = gsub("ADT-", "", data_names)) %>%
    tidyr::separate(TEMP, into = c("TEMP", "Clone")) %>%
    dplyr::mutate(TEMP = toupper(gsub(" ", "", TEMP)),
                  TEMP = case_when(TEMP == "HLA" ~ "HLA-DR",
                                   TEMP == "PD1L1" ~ "CD274",
                                   TRUE ~ TEMP))

witkowski_2020 <- citeseq %>%
    dplyr::filter(Study == nm) %>%
    dplyr::mutate(TEMP = toupper(gsub(" ", "", Antigen)))

setdiff(data_names$TEMP, witkowski_2020$TEMP)
setdiff(witkowski_2020$TEMP, data_names$TEMP)

witkowski_2020 <- witkowski_2020 %>%
    full_join(data_names %>% select(-Clone), by = "TEMP") %>%
    tidyr::fill(Study) %>%
    dplyr::mutate(Antigen = dplyr::coalesce(Antigen, TEMP)) %>%
    dplyr::select(-TEMP)

# Wu_2021 ----

# Data is missing some of the reported antibodies
# By private communication, the authors told us that either -
#  - the antibody was included but had no counts
#  - the antibody was not included but was present through contamination

# It is not clear from the data whether the isotype controls are mouse or rat
# as some are missing

wu_patch <- data.frame(Data_Name = c("MHCII-AATAGCGAGCAAGTA",
                                     "MHCI-TATGCGAGGCTTATC",
                                     "IgG2b-ATATGTATCACGCGA",
                                     "FCER1A-CTCGTTTCCGTATCG",
                                     "CD40LG-GCTAGATAGATGCAA",
                                     "BT3_1-TACTTCTGGAGCCAG",
                                     "CR2-AACCTAGTAGTTCGG",
                                     "TCR_V24-J18-AACTTCTGTGGTAGC",
                                     "CD8a-GCTGCGCTTTCCATT"),
                         TEMP = c("HLADR", "HLAABC", "MOUSE IGG2BISO",
                                  "FCERIA", "CD40L", "CD277", "CD21",
                                  "TCRVA24JA18", "CD8A"))

nm <- "Wu_2021_b"
print(nm)

#data_names <- unique(unlist(lapply(ab_fnames[["Wu2021b"]], read_rds)))
data_names <- protein_names[[acc_to_study[[nm]]]]

data_names <- data.frame(Data_Name = data_names) %>%
    dplyr::filter(! Data_Name == "unmapped") %>%
    tidyr::separate(Data_Name, into = c("TEMP", "Barcode_Sequence"),
                    remove = FALSE, sep = "-(?=[ACTG]+)") %>%
    dplyr::mutate(TEMP = toupper(gsub("[-_]", "", TEMP))) %>%
    dplyr::rows_update(wu_patch, by = "Data_Name")

wu_2021 <- citeseq %>%
    dplyr::filter(Study == "Wu_2021") %>%
    dplyr::mutate(TEMP = toupper(gsub("\\.", "", Antigen)))

# Just the isotype controls unmatched in data
setdiff(wu_2021$TEMP, data_names$TEMP)
setdiff(data_names$TEMP, wu_2021$TEMP)

wu_2021 <- wu_2021 %>%
    dplyr::left_join(data_names %>% select(-Barcode_Sequence), by = "TEMP") %>%
    dplyr::rows_patch(data_names, unmatched = "ignore") %>%
    dplyr::select(-TEMP)

# 10x ----

print("10x")

tenx_data_names <- protein_names[grepl("10x",names(protein_names))]
#tenx_data_names <- lapply(tenx_fnames, read_rds)
tenx_data_names <- data.frame(Data_Name = unlist(tenx_data_names),
                              Study = rep(names(tenx_data_names),
                                          lengths(tenx_data_names))) %>%
    dplyr::mutate(TEMP = gsub("-Total(Seq)?[ABCD]$", "", Data_Name))

tenx_clones <- citeseq %>%
    dplyr::filter(grepl("10x", Study)) %>%
    dplyr::mutate(TEMP = paste(Antigen, Clone, sep = "-"),
                  TEMP = gsub("_", "-", TEMP)) %>%
    # Some data names include Clone, first add these
    dplyr::left_join(tenx_data_names, by = c("Study", "TEMP"))

# Find the unmatched Data_Names
tenx_data_names <- tenx_data_names %>%
    dplyr::anti_join(tenx_clones) %>%
    dplyr::mutate(TEMP = gsub("-control", "", TEMP))
                  #TEMP = ifelse(Study %in% c("10x19Nov2018", "10x19Nov2018-2") &
                  #                  grepl("control", Data_Name),
                  #              gsub("-control", "", TEMP), TEMP))

tenx_clones <- tenx_clones %>%
    dplyr::mutate(TEMP = Antigen) %>%
    dplyr::rows_patch(tenx_data_names, by = c("TEMP", "Study"))


# Check that all data names are correct:
anti_join(tenx_data_names, tenx_clones)

# Combine data names, patch into citeseq data -----

combined <- dplyr::bind_rows(arunachalam_2020, hao_2021,
    buus_2021, cadot_2020, chung_2021, fernandez_2020, fidanza_2020,
    frangieh_2021, granja_2019, holmes_2020, kaufmann_2021, kotliarov_2020,
    krebs_2020, lawlor_2021, leader_2021, lecoz_2021, lee_2021, liu_2021,
    mair_2020, mimitou_2019, mimitou_2021, nathan_2021,
    papalexi_2021, pei_2020, poch_2021, pomboAntunes_2021, pont_2020, qi_2020,
    rinconArevalo_2021, shangguan_2021, stephenson_2021, stoeckius_2017,
    stoeckius_2018,
    stuart_2019, su_2020, trzupek_2020, trzupek_2021, valenzi_2019, wang_2020,
    vanuytsel_2020, witkowski_2020, wu_2021, tenx_clones) %>%
    dplyr::mutate(Antigen = coalesce(Antigen, Data_Name)) %>%
    dplyr::select(-ID) %>%
    AbNames::addID()


combined_add <- dplyr::anti_join(citeseq, combined, by = c("ID"))

citeseq <- combined %>%
    dplyr::bind_rows(combined_add) %>%
    dplyr::select(-ID, -TEMP) %>%
    dplyr::relocate(ID, Antigen, Data_Name) %>%
    AbNames::addID()

#readr::write_delim(citeseq,
#                   file = "~/Analyses/AbNames/inst/extdata/citeseq.csv",
#                   delim = ",")

#readr::write_delim(citeseq, file = "ADT_clones/merged_adt_clones.tsv")


#"Arunachalam_2020" "Stuart_2019"
# "Granja_2019"      "Hao_2021" "Stoeckius_2018"

# To do - add a column indicating if in data not in clone table?
# Check if antigens missing from combined table are correct
