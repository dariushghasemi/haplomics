#!/usr/bin/Rscript


# The script was created on 22 December 2023
# to automate the process of preparing data
# for haplotype reconstruction analysis on
# all of the 11 replicated kidney loci.

#----------#

is.installed <- function(package_name){
    is.element(package_name, installed.packages()[,1])
}

# check if the package is installed
if (!is.installed("config")){
   install.packages("config");
}

#library(config)
library(tidyverse)

config <- config::get(file = "../config/configuration.yml")

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
genotype_file <- args[1]
ophen <- args[2]
ometa <- args[3]
oprot <- args[4]

# data files imported from config file
traits_blood   <- config$path_phen
traits_metabol <- config$path_meta
traits_protein <- config$path_prot
principal_comp <- config$path_comp


# taking the locus name
locus_name  <- gsub("_dosage.txt", "", basename(genotype_file))
locus_name

#----------#
# directories

#out.phen.rds <- paste0(out.dir, locus_name, "_haplotypes_data_phen.RDS")
#out.meta.rds <- paste0(out.dir, locus_name, "_haplotypes_data_meta.RDS")
#out.prot.rds <- paste0(out.dir, locus_name, "_haplotypes_data_prot.RDS")


#-----------------------------------------------------#
#-------               Functions             ---------
#-----------------------------------------------------#

# inverse-normal transformation
do_INT = function(values) {
  qnorm((rank(values, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(values)))
}

# median imputation
median_imput = function(x) ifelse(is.na(x), median(x, na.rm = T), x)

# quantitative traits in the CHRIS study for PheWAS
phenotypes <- c(
    "Height","Weight","BMI","Body_Fat","Visceral_Fat",
    "SBP", "DBP", "Pulse_Rate", "HbA1c",
    "SAlb", "SCr", "UAlb", "UCr", "UACR", "eGFRw",
    "PTT", "INR_PT", "APTT_ratio", "APTT", "Fibrinogen", "AT", "BGlucose",
    "Urate", "AST_GOT", "ALT_GPT", "GGT", "ALP", "TB", "DB", "Lipase", "TC",
    "HDL", "LDL", "TG", "Sodium", "Potassium", "Chlorine", #"Calcium_mg",
    "Calcium", "Phosphorus", "Magnesium", "Iron", "Ferritin", 
    "Transferrin", "TIBC" , "TS", "Homocyst", "CRP", "TSH", #"FT3", "FT4", 
    "Cortisol", "WBC","RBC","HGB","HCT", "MCV","MCH","MCHC","RDW","PLT","MPV",
    "Neutrophils","Lymphocytes","Monocytes","Eosinophils","Basophils",
    "AntiTPO","Urine_pH","UGlucose","UProteins","UHGB",
    "S2UCr", "S2UAlb", "UG2Cr", "FE_Glu", "FE_Alb", "FE_HGB"
  )


#-----------------------------------------------------#
#-------              Import data            ---------
#-----------------------------------------------------#

cat("\nImport data...\n")

# Genotype data
genome <- read.delim(
  genotype_file,
  col.names  = c("AID", "chromosome", "position", "MARKER_ID", "REF", "ALT", "AF", "Dosage"),
  colClasses = c(rep("character", 6), rep("numeric", 2)),
  #nrows = 125000
)

# CHRIS metabolites; p = 175; n = 7,252
# defining AID column as character to preserve 
# the leading zero in participants identifier
chrisMass <- read.csv(traits_metabol, colClass = c(AID = "character")) %>%
  # Excluding with missing rate > 10.2%
  select(- c(pc_aa_c30_2, met_so, pc_ae_c38_1))

# Saving metabolites names for PheWS
metabolites <- chrisMass %>% select(- AID) %>% colnames()

# CHRIS proteins concentrations; p = 148; n = 4,087 
chrisProt <- read.csv(traits_protein, colClass = c(AID = "character"))

# Saving proteins names for PheWS
proteins <- chrisProt %>% select(- AID) %>% colnames()

# Principle components
prc_comp <- read.delim(principal_comp, sep = "\t", stringsAsFactors = FALSE) %>% rename(AID = X.IND_ID)

# Phenotype data 0 previously used pheno: 06-Nov-2022_chris4HaploReg.txt
chris <- read.csv(traits_blood, header =TRUE, stringsAsFactors = FALSE)

#----------#
# making new phenotypes
chris <- chris %>%
  as_tibble() %>%
  mutate(
    # Serum-to-Urine-Creatinine ratio
    S2UCr  = SCr / UCr,
    # Serum-to-Urine-Albumin ratio
    S2UAlb = SAlb / UAlb,
    # Urinary-Glucose-to-Creatinine ratio
    UG2Cr  = UGlucose / UCr,
    # fractional excretion of Glucose levels
    FE_Glu = (UGlucose * SCr) / (BGlucose * UCr),
    # fractional excretion of Albumin
    FE_Alb = (UAlb * SCr) / (SAlb * UCr),
    # fractional excretion of Humoglobin
    FE_HGB = (UHGB * SCr) / (HGB * UCr)
)

# Saving traits names for PheWAS
phenome <- chris %>% select(all_of(phenotypes)) %>% colnames()

#-----------------------------------------------------#
#-------              Dictionary             ---------
#-----------------------------------------------------#


#-----------------------------------------------------#
#-------               Merge data            ---------
#-----------------------------------------------------#

cat("\nMerge genotypes...\n")

# Restructuring vcf file to wide format and merged with PCs
to_merge_genome <- genome %>%
  filter(AF > 0.00001) %>%
  mutate(SNPid = str_c("chr", chromosome, ":", position)) %>%
  distinct(AID, SNPid, .keep_all = TRUE) %>%
  pivot_wider(id_cols = AID, names_from = SNPid, values_from = Dosage) %>%
  inner_join(by = "AID", prc_comp) %>%
  inner_join(by = "AID", chris %>% select(AID, Sex, Age)) %>%
  mutate(
    across(Age, as.numeric),
    across(Sex, as.factor)
  )

str(to_merge_genome)
dim(to_merge_genome)
str(to_merge_genome[c("AID", "Age", "Sex")])

#----------#
cat("\nMerged clinical traits...\n")

# dataset containing genotypes and clinical traits
merged_phen <- to_merge_genome %>%
  inner_join(by = "AID", chris %>% select(AID, all_of(phenome))) %>%
  mutate(
    across(any_of(phenome), as.numeric),
    across(any_of(phenome), median_imput)
    ) #%>% slice_head(n = 6000) 

str(merged_phen %>% select(- starts_with("chr")))

cat("\nMerged serum metabolites...\n")

# dataset containing genotypes and serum metabolites
merged_meta <- to_merge_genome %>%
  inner_join(by = "AID", chrisMass) %>%
  mutate(across(any_of(metabolites), median_imput))

str(merged_meta %>% select(- starts_with("chr")))

cat("\nMerged plasma proteins...\n")

# dataset containing genotypes and plasma proteins
merged_prot <- to_merge_genome %>%
  inner_join(by = "AID", chrisProt) %>%
  mutate(across(any_of(proteins), median_imput))

str(merged_prot %>% select(- starts_with("chr")))

#----------#
# save the merged data
write.csv(merged_phen, file = ophen, quote = F, row.names = F)
write.csv(merged_meta, file = ometa, quote = F, row.names = F)
write.csv(merged_prot, file = oprot, quote = F, row.names = F)

#saveRDS(merged_data, file = output.rds)

#-----------------------------------------------#
# Printing number of SNPs and traits
message.data <- paste(
  "Store", 
  locus_name, 
  "haplotypes data including dosage levels of",
  ncol(to_merge_genome %>% select(starts_with("chr"))),
  "variants and",
  length(intersect(phenome, colnames(merged_phen))),
  "clinical traits\nfor ", 
  nrow(merged_phen),
  "individuals and",
  length(intersect(metabolites, colnames(merged_meta))),
  "serum metabolites for", 
  nrow(merged_meta),
  "individuals and",
  length(intersect(proteins, colnames(merged_prot))),
  "plasma proteins (n =",
  nrow(merged_prot),
  ") individuals on: \n",
  ophen)

cat("\n", message.data, "\n")

#----------#
# print time and date
Sys.time()

#sbatch --wrap 'Rscript 03-1_haplotypes_data.R data/dosage/PDILT_dosage.txt' -c 2 --mem-per-cpu=16GB -J "03-1_PDILT.R"