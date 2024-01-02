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
pheno.geno <- args[1]

# taking the locus name
locus_name  <- gsub("_haplotypes_data.csv", "", basename(pheno.geno))
locus_name

#----------#
# directories

base.dir   <- "/home/dghasemisemeskandeh/projects/haploAnalysis/output/result_associations/"
data.dir   <- paste0(base.dir, "/data/pheno/")
out.dir    <- paste0(base.dir, "/output/result_association/")
#pheno.geno <- paste0(data.dir, locus_name, "_haplotypes_data.RDS")
output.rds <- paste0(base.dir, locus_name, "_haplotypes_association_min1.RDS")


#-----------------------------------------------------#
#-------               Functions             ---------
#-----------------------------------------------------#

# inverse-normal transformation
do_INT = function(values) {
  qnorm((rank(values, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(values)))
}


#-----------------------------------------------------#
#-------              Import data            ---------
#-----------------------------------------------------#

cat("\nImport data...\n")

# Phenotype data merged with genotypes
merged_data <- read.csv(pheno.geno, header = TRUE, stringsAsFactors = FALSE)


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
haplo_dataset  <- data.frame(haplo_genotype, merged_data %>% select(APTT,AST_GOT,Cortisol,DBP,SBP,Pulse_Rate,eGFRw, Age, Sex, PC1:PC10)) #(-AID, -starts_with("chr")))

#merged_data %>% select(TSH, eGFRw, HDL, APTT, BMI, Age, Sex, PC1:PC10))
#(ALT_GPT,FE_Alb,Iron,TS, eGFRw, Age, Sex, PC1:PC10)) #
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
  
  # Set parameters to control EM algorithm
  em_ctrl <- haplo.em.control(
    n.try = 2,
    iseed = common_iseed,
    insert.batch.size = 2, # keep it =2 to equalize n.of haplo for all traits
    max.haps.limit = 4e6,
    min.posterior = 1e-1  # increase the prob of trimming off rare haplo at each insertion step (raise of prob will end up qith fewer haplos)
    )
  
  # fiting the model
  model_fit <- haplo.glm(
    my_formula,
    family = gaussian,
    data   = df,
    na.action = "na.geno.keep",
    locus.label = colnames(loci),
    x = TRUE,
    control = haplo.glm.control(
      haplo.freq.min = .01, # if drop it to 0.01, then it unequlizes no. haplotypes
      em.c = em_ctrl
      )
    )
  
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
saveRDS(results_shrinked, output.rds)

#----------#
# print time and date
Sys.time()

#sbatch --wrap 'Rscript 03-2_haplotypes_building.R data/pheno/IGF1R_haplotypes_data.csv' -c 2 --mem-per-cpu=32GB -J "03-2_IGF1R.R"