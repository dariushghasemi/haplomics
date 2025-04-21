
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--dosage", default=NULL, help="Path and filename of master coloc table produced by individual traits pre-processing"),
  make_option("--phenotype", default=NULL, help="Path to phenotypes file"),
  make_option("--covariate", default=NULL, help="Path to covariates file"),
  make_option("--min_ac", default=1, help="Minimum allele count threshold to remove mono-allelic variants with mac below it"),
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

#----------#
# number of individuals in genotype file
n_snps <- length(unique(paste0(genome$CHR,":",genome$POS)))
n_sample <- length(unique(genome$IID))
min_ac <- opt$min_ac
max_ac <- n_sample - opt$min_ac


cat("\nReshape genotypes to wide...\n")

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
cat(
  "\nOf",
  n_snps,
  "variants,",
  ncol(genome_wide) - 1,
  "left with minor allele count >",
  opt$min_ac,
  "for building haplotypes.\n"
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


