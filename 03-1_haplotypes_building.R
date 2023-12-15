#!/usr/bin/Rscript

# The script was created on 12 December 2023
# to automate haplotype reconstruction analysis
# for all of the 11 replicated kidney loci.

#----------#
library(tidyverse)
library(haplo.stats)

#----------#

is.installed <- function(package_name){
    is.element(package_name, installed.packages()[,1])
}

# check if package "Haplo.stats" is installed
if (!is.installed("haplo.stats")){
   install.packages("haplo.stats");
}

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

base.dir <- "/home/dghasemisemeskandeh/projects/haploAnalysis"
out.dir  <- paste0(base.dir, "/output/result_associations")
traits_blood   <- "/home/dghasemisemeskandeh/projects/HaploReg/data/chris_q-norm.csv"
principal_comp <- "/home/dghasemisemeskandeh/projects/HaploReg/data/CHRIS13K.GT.evecs"
output.full    <- paste0(out.dir, "/", today.date, "_", locus_name, "_association_results_full1.RDS")
output.short   <- paste0(out.dir, "/", today.date, "_", locus_name, "_association_results_short1.RDS")


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
    "Transferrin", "TIBC" , "TS", "Homocyst", "CRP", "TSH", "FT3", "FT4", 
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
  pivot_wider(id_cols = AID, names_from = SNPid, values_from = Dosage) %>%
  inner_join(by = "AID", chris %>% select(
    AID,
    Sex,
    Age,
    eGFRw,
    all_of(phenome[c(3,6,7,9)]),
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


#-----------------------------------------------------#
#-------            Haplotype data           ---------
#-----------------------------------------------------#

cat("\nBuilding data...\n")

# S1: selecting the variants
loci <- merged_data %>% select(starts_with("chr"))

#----------#
# Number of variants
cat("\n", ncol(loci), "SNPs used for building haplotypes at", locus_name, "locus.", "\n")
#----------#

# S2: convert the dosage to integer, then to major(1/2) or minor(1/2) alleles
genome_binary  <- geno1to2(round(loci, 0), locus.label = colnames(loci))

# S3: setup genotype data
haplo_genotype <- setupGeno(genome_binary, miss.val = c(0, NA), locus.label = colnames(loci))

# S4: GLM data (merging phenotype and genotype data)
haplo_dataset  <- data.frame(haplo_genotype,
                             merged_data %>% select(-AID, -starts_with("chr")))

#----------#
dim(haplo_dataset)
str(haplo_dataset)

#-----------------------------------------------------#
#-------       Haplotype reconstruction      ---------
#-----------------------------------------------------#

# Fitting regression model
# glm fit with haplotypes, additive gender covariate on gaussian response
#----------#

# Defining Haplo.GLM model for iteration via map function
hap_model <- function(df){
  
  # making PCs vector to adjust
  PCs <- paste0("PC", 1:10, collapse = " + ")
  
  # defining model formula
  my_formula <- paste("do_INT(trait) ~ haplo_genotype + Sex + Age +", PCs)

  # Set a common random number seed for both models
  common_iseed <- 777
  
  # Set haplo.em.control parameters
  em_ctrl <- haplo.em.control(
    n.try = 2,
    iseed = common_iseed,
    insert.batch.size = 2,
    max.haps.limit = 4e6,
    min.posterior = 1e-5)
  
  # fiting the model
  model_fit <- haplo.glm(
    my_formula,
    family = gaussian,
    data   = df,
    na.action = "na.geno.keep",
    locus.label = colnames(loci),
    x = TRUE,
    control = haplo.glm.control(haplo.freq.min = .02,
                                em.c = em_ctrl))
  
  return(model_fit)
}

#----------#
# Retrieving haplotypes from fitted models and their frequencies
haplo_extract <- function(model) {
  
  haplo_set <- 
    summary(model)$haplotypes %>%
    rownames_to_column(var = "Haplotype") %>%
    mutate(Haplotype = str_replace(
      Haplotype,
      "(?<=\\.)\\d{1,2}(?!\\d)",
      sprintf("%03d", as.numeric(str_extract(Haplotype, "(?<=\\.)\\d{1,2}(?!\\d)")))),
      Haplotype = str_replace_all(Haplotype,
                                  c("haplo_genotype." = "H",
                                    "haplo.base"      = "Ref.")))
  
  return(haplo_set)
}


#-----------------------------------------------------#
#-------        Haplotype association        ---------
#-----------------------------------------------------#

cat("\nBuilding model...\n")

# Iterating the model on the traits
results <- haplo_dataset %>%
  #select(- any_of(quantVars[c(1:13, 15:38, 40:76)])) %>%
  pivot_longer(cols      = - c(haplo_genotype, Age, Sex, PC1:PC10),
               names_to  = "trait_name",
               values_to = "trait") %>%
  group_by(trait_name) %>%
  nest() %>%
  mutate(model     = data  %>% map(hap_model),
         haplotype = model %>% map(haplo_extract),
         glance    = model %>% map(broom::glance),
         tidy      = model %>% map(broom::tidy))


results %>% unnest(tidy)

#----------#
# saving full model results

cat("\nSaving results...\n")

# Drop unnecessary results
results_shrinked <- results %>% select(trait_name, haplotype, tidy)

# saving the results
saveRDS(results, output.full)
saveRDS(results_shrinked, output.short)

#----------#
# print time and date
Sys.time()

#sbatch --wrap 'Rscript 03-1_haplotypes_building.R genotype/PDILT_dosage.txt' -c 2 --mem-per-cpu=16GB -J "03-1_PDILT.R"