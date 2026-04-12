
#-----------------------------------------------------#
#-------                Logging              ---------
#-----------------------------------------------------#

# Print time & message each time
log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(timestamp, " - ", ...)
}

#-----------------------------------------------------#
#------              Load Libraries            -------
#-----------------------------------------------------#

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
})

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--dosage", type="character", help="Dosage file"),
  make_option("--phenotype", type="character", help="Phenotypes file"),
  make_option("--covariate", type="character", help="Covariates file (optional)"),
  make_option("--min_ac", type="numeric", default=NULL, help="Minimum allele count threshold to remove mono-allelic variants with mac below it"),
  make_option("--summary", type="character", help="Summary filename"),
  make_option("--output", type="character", help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


log_message("Parsed parameters: ", opt)


#-----------------------------------------------------#
#-------              Import data            ---------
#-----------------------------------------------------#

log_message("Importing data...")


# phenotypes
pheno_file <- data.table::fread(opt$phenotype, colClass = c(IID = "character"))

# list of traits in phenotype file
phenome <- colnames(pheno_file %>% dplyr::select(- IID))


# column names in dosage file
dosage_cols <- c("IID", "CHR", "POS", "MARKER_ID", "REF", "ALT", "AF", "Dosage")
dosage_type <- c(rep("character", 6), rep("numeric", 2))


# Genotype data
genome <- read.delim(
  opt$dosage,
  col.names  = dosage_cols,
  colClasses = dosage_type,
  na.strings = "."
)

#----------#
# number of individuals in genotype file
n_snps <- length(unique(paste0(genome$CHR,":",genome$POS)))
n_sample <- length(unique(genome$IID))
min_ac <- opt$min_ac
max_ac <- n_sample - opt$min_ac

rep_l1 <- paste(
  "Minimum allele count (AC) is set to", opt$min_ac,
  ". Given", n_sample,
  "samples, variants with AC below", min_ac,
  "and above", max_ac,
  "are removed."
  )

log_message("✅ Done! ", rep_l1)

#----------#
log_message("Reshape genotypes to wide table...")

# Restructuring vcf file to wide format and merged with PCs
genome_wide <- genome %>%
  dplyr::mutate(
    SNPid = paste0("chr", CHR, "_", POS, "_", REF, "_", ALT),
    MAF = ifelse(AF > 0.5, 1 - AF, AF),
    AC = AF * n_sample
    ) %>%
  dplyr::filter(
    #AF > 0.001 & AF < 0.999  # filter out mono-allelic variants respect to H-W equilibrium
    AC > min_ac & AC < max_ac
    ) %>%
  pivot_wider(
    id_cols = IID, 
    names_from = SNPid, 
    values_from = Dosage
    )

# show how nany variants filtered out for AC limit
rep_l2 <- paste(
  "Of",
  n_snps,
  "variants,",
  ncol(genome_wide) - 1,
  "left with minor allele count >",
  opt$min_ac,
  "for building haplotypes."
  )

log_message("✅ Done! ", rep_l2)


#-----------------------------------------------------#
#-------               Merge data            ---------
#-----------------------------------------------------#

log_message("Merging data...")

if(!is.null(opt$covariate) && opt$covariate != "" && opt$covariate != "None"){
  
  # Principle components and covariates
  covar_file <- data.table::fread(opt$covariate, colClass = c(IID = "character"))
  
  merged_file <- genome_wide %>%
    inner_join(by = "IID", pheno_file) %>%
    inner_join(by = "IID", covar_file) %>%
    dplyr::mutate(
      across(any_of(phenome), as.numeric)
      )
  
  rep_l3 <- paste(
    "Merged data comprised",
    ncol(genome_wide[-1]), "variants,", 
    ncol(pheno_file[-1]), "traits, and", 
    ncol(covar_file[-1]), "covariates."
    )
  
  log_message("✅ Done! ", rep_l3)
  
  } else{
  # dataset containing genotypes and clinical traits
  merged_file <- genome_wide %>%
    inner_join(by = "IID", pheno_file) %>%
    dplyr::mutate(
      across(any_of(phenome), as.numeric)
    )
  
  rep_l3 <- paste(
    "Merged data comprised",
    ncol(genome_wide[-1]), "variants,", 
    ncol(pheno_file[-1]), "traits. No covariate file provided."
  )
  
  log_message("✅ Done! ", rep_l3)
  }

#----------#
log_message("Saving the report...")

# Combine and save report to file
full_report <- paste(rep_l1, rep_l2, rep_l3, sep = " ")

# Save report to file
writeLines(full_report, con = opt$summary)

log_message("✅ Done! Find it here: ", opt$summary)

#----------#
log_message("Saving merged data...")

# save the merged data
saveRDS(merged_file, file = opt$output)

log_message("✅ Done! Find it here: ", opt$output)
