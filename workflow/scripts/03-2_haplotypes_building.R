# load libraries
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(haplo.stats, quietly = TRUE))
suppressMessages(library(optparse, quietly = TRUE))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--data", default=NULL, help="Merged genotype, phenotypes, and covariates file"),
  make_option("--covariate", default=NULL, help="Covariates file (optional)"),
  make_option("--min_freq", default=0.01, help="Minimum frequency for rare haplotypes"),
  make_option("--max_haps", default=4e6, help="Maximum number of haplotypes"),
  make_option("--min_pp", default=1e-5, help="Minimum posterior probability"),
  make_option("--n_batch", default=2, help="Number of batches"),
  make_option("--n_try", default=2, help="Number of attempts"),
  make_option("--output", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#----------#
# print time and date
Sys.time()

# date
today.date <- format(Sys.Date(), "%d-%b-%y")


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
merged_data <- readRDS(opt$data)


# checking whether the covariate file was provided and not empty
if(!is.null(opt$covariate) && opt$covariate != "" && opt$covariate != "None"){ 
  
  # covariate file
  covar_data <- data.table::fread(opt$covariate)
  
  # Extract covariate names excluding ids
  covariates <- colnames(covar_data[, !"IID", with = FALSE])
  
  # Build covariate terms to adjust
  covar_term <- paste0(covariates, collapse = " + ")
  
  # defining model formula with covariates
  model_formula <- paste("do_INT(trait) ~ haplo_genotype +", covar_term)
  
} else{
  
  # defining model formula without covariates
  model_formula <- paste("do_INT(trait) ~ haplo_genotype")
  covariates <- NULL

}



#-----------------------------------------------------#
#-------            Haplotype data           ---------
#-----------------------------------------------------#

cat("\nBuilding data...\n")

# S1: selecting the variants
loci <- merged_data %>% dplyr::select(matches("chr[0-9]{,2}_[0-9]+_[ATCG]+_[ATCG]+"))

# S2: convert the dosage to integer, then to major(1/2) or minor(1/2) alleles
genome_binary  <- haplo.stats::geno1to2(round(loci, 0), locus.label = colnames(loci))

# S3: setup genotype data
haplo_genotype <- haplo.stats::setupGeno(genome_binary, miss.val = c(0, NA), locus.label = colnames(loci))

# S4: GLM data (merging phenotype and genotype data)
haplo_dataset  <- data.frame(haplo_genotype, merged_data %>% dplyr::select(-IID, -matches("chr[0-9]{,2}_[0-9]+_[ATCG]+_[ATCG]+")))

#----------#
# Number of variants
cat("\n", ncol(loci), "SNPs used for building haplotypes at this locus.", "\n")


#-----------------------------------------------------#
#-------       Haplotype reconstruction      ---------
#-----------------------------------------------------#

# Fitting regression model: additive haplotypes and covariate on gaussian response
# Defining Haplo.GLM model for iteration via map function
# recall parameters of haplo.GLM model from opt arguments
hap_model <- function(
    df,
    n_try = opt$n_try,
    n_batch  = opt$n_batch,
    max_haps = opt$max_haps,
    min_pp   = opt$min_pp,
    min_freq = opt$min_freq,
    my_formula = model_formula,
    common_iseed = 777   # Set a common random number seed for both models
    ){
  
  
  # Set parameters to control EM algorithm
  em_ctrl <- haplo.stats::haplo.em.control(
    n.try = n_try,            # to have equal no. of haplotypes, phenotype must be imputed! 
    iseed = common_iseed,
    insert.batch.size = n_batch, # keep it =2 to equalize n.of haplo for all traits
    max.haps.limit = max_haps,
    min.posterior = min_pp  # increase the prob of trimming off rare haplo at each insertion step (raise of prob will end up qith fewer haplos)
    )
  
  # fiting the model
  model_fit <- haplo.stats::haplo.glm(
    formula = my_formula,
    family = gaussian,
    data   = df,
    na.action = "na.geno.keep",
    locus.label = colnames(loci),
    x = TRUE,
    control = haplo.glm.control(
      haplo.freq.min = min_freq, # if drop it to 0.01, then it unequlizes no. haplotypes
      em.c = em_ctrl
      )
    )
  
  return(model_fit)
}

#----------#
# Retrieving haplotypes from fitted models and their frequencies
haplo_extract <- function(model) {
  
  haplo_set <- summary(model)$haplotypes %>%
    tibble::rownames_to_column(var = "Haplotype") %>%
    dplyr::mutate(
      Haplotype = str_replace(
        Haplotype,
        "(?<=\\.)\\d{1,2}(?!\\d)",
        sprintf("%03d", as.numeric(str_extract(Haplotype, "(?<=\\.)\\d{1,2}(?!\\d)")))
        ),
      Haplotype = str_replace_all(
        Haplotype, c(
          "haplo_genotype." = "H", 
          "haplo.base"      = "Ref."
          )
        )
      )
  
  return(haplo_set)
}


#-----------------------------------------------------#
#-------        Haplotype association        ---------
#-----------------------------------------------------#

cat("\nBuilding model...\n")

# Iterating the model on the traits
results <- haplo_dataset %>%
  pivot_longer(
    cols      = - c(haplo_genotype, all_of(covariates)),
    names_to  = "trait_name",
    values_to = "trait"
    ) %>%
  group_by(trait_name) %>%
  nest() %>%
  dplyr::mutate(
    model     = data  %>% map(hap_model),
    haplotype = model %>% map(haplo_extract),
    glance    = model %>% map(broom::glance),
    tidy      = model %>% map(broom::tidy)
    ) %>%
  ungroup()


results %>% unnest(tidy)

#----------#
# saving full model results

cat("\nSaving results...\n")

# Drop unnecessary results
results_shrinked <- results %>% dplyr::select(- c(data, model, glance))

# saving the results
saveRDS(results_shrinked, opt$output)

#----------#
# print time and date
Sys.time()
