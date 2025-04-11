
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--dosage", default=NULL, help="Path and filename of master coloc table produced by individual traits pre-processing"),
  make_option("--phenotype", default=NULL, help="Path to phenotypes file"),
  make_option("--covariate", default=NULL, help="Path to covariates file"),
  make_option("--mac", default=20, help="Minor allele count threshold to filter out variants below it"),
  make_option("--output", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#-----------------------------------------------------#
#-------               Functions             ---------
#-----------------------------------------------------#

# inverse-normal transformation
do_INT = function(values) {
  qnorm((rank(values, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(values)))
}

# median imputation
median_imput = function(x) ifelse(is.na(x), median(x, na.rm = T), x)



#-----------------------------------------------------#
#-------              Import data            ---------
#-----------------------------------------------------#

cat("\nImport data...\n")


# phenotypes
pheno_file <- data.table::fread(opt$phenotype) #colClass = c(IID = "character")

# list of traits in phenotype file
phenome <- colnames(pheno_file %>% dplyr::select(- IID))

# Principle components and covariates
covar_file <- data.table::fread(opt$covariate)

#----------#

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


cat("\nReshape genotypes to wide...\n")

# Restructuring vcf file to wide format and merged with PCs
genome_wide <- genome %>%
  #filter(AF > 0.00001) %>%
  dplyr::mutate(
    SNPid = paste("chr", CHR, POS, REF, ALT, sep = "_")
    ) %>%
  #distinct(AID, SNPid, .keep_all = TRUE) %>%
  pivot_wider(
    id_cols = IID, 
    names_from = SNPid, 
    values_from = Dosage
    )


#-----------------------------------------------------#
#-------               Merge data            ---------
#-----------------------------------------------------#

cat("\nMerge with traits...\n")

if(!file.exists(opt$covariate)){
  
  # dataset containing genotypes and clinical traits
  merged_file <- genome_wide %>%
    inner_join(by = "IID", pheno_file) %>%
    dplyr::mutate(
      across(any_of(phenome), as.numeric),
      across(any_of(phenome), median_imput)
    ) #%>% slice_head(n = 6000)
  
  } else{
    
    merged_file <- genome_wide %>%
      inner_join(by = "IID", pheno_file) %>%
      inner_join(by = "IID", covar_file) %>%
      dplyr::mutate(
        across(any_of(phenome), as.numeric),
        across(any_of(phenome), median_imput),
        across(Age, as.numeric),
        across(Sex, as.factor)
      )
    }


#----------#
cat("\nSave output datasets...\n")

# save the merged data
saveRDS(merged_file, file = opt$output)


