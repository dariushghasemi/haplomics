#!/usr/bin/Rscript


# The script was created on 22 December 2023
# to automate the process of preparing data
# for haplotype reconstruction analysis on
# all of the 11 replicated kidney loci.

#----------#

library(tidyverse)

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")

#----------#
# taking variants file as input
args <- commandArgs(trailingOnly = TRUE)
genotype_file <- args[1]

# taking the locus name
locus_name  <- gsub("_dosage.txt", "", basename(genotype_file))
locus_name

#----------#
# directories

base.dir <- "/home/dghasemisemeskandeh/projects"
out.dir  <- paste0(base.dir, "/haploAnalysis/data/pheno/")
data.dir <- paste0(base.dir, "/HaploReg/data")
traits_blood   <- paste0(data.dir, "/chris_q-norm.csv")
principal_comp <- paste0(data.dir, "/CHRIS13K.GT.evecs")
output.csv <- paste0(out.dir, locus_name, "_haplotypes_data.csv")
output.rds <- paste0(out.dir, locus_name, "_haplotypes_data.RDS")


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
    "SAlb", "SCr", "UAlb", "UCr", "UACR",
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
#-------               Merge data            ---------
#-----------------------------------------------------#

cat("\nMerge data...\n")

# Restructuring vcf file to wide format
merged_data <- genome %>%
  filter(AF > 0.00001) %>%
  mutate(SNPid = str_c("chr", chromosome, ":", position)) %>%
  distinct(AID, SNPid, .keep_all = TRUE) %>%
  pivot_wider(id_cols = AID, names_from = SNPid, values_from = Dosage) %>%
  inner_join(by = "AID", chris %>% select(
    AID,
    Sex,
    Age,
    eGFRw,
    all_of(phenome),
    #-FT3,
    #-FT4
    )
    ) %>%
  mutate(
    across(c("Age", any_of(phenome)), as.numeric),
    across(c("Sex"), as.factor)
    ) %>%
  #slice_head(n = 6000) %>%
  inner_join(prc_comp, by = "AID")

str(merged_data)
#----------#
# save the merged data
write.csv(merged_data, file = output.csv, quote = F, row.names = F)
saveRDS(merged_data, file = output.rds)

#-----------------------------------------------#
# Printing number of SNPs and traits
message.data <- paste(
  "Store", 
  locus_name, 
  "haplotypes data including dosage levels of",
  ncol(merged_data %>% select(starts_with("chr"))),
  "variants and",
  length(intersect(phenome, colnames(merged_data))),
  "traits on: \n",
  output.rds)

cat("\n", message.data, "\n")

#----------#
# print time and date
Sys.time()

#sbatch --wrap 'Rscript 03-1_haplotypes_data.R data/dosage/PDILT_dosage.txt' -c 2 --mem-per-cpu=16GB -J "03-1_PDILT.R"