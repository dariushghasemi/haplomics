
#-----------------------------------------------------#
# ------             Load libraries             ------
#-----------------------------------------------------#

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(haplo.stats)
  library(optparse)
})

#----------#
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
# params
i_seed   = 777
n_try    = opt$n_try
n_batch  = opt$n_batch
max_haps = opt$max_haps
min_pp   = opt$min_pp
min_freq = opt$min_freq
ofile = opt$output
cfile = opt$covariate


#-----------------------------------------------------#
#-------                Logging              ---------
#-----------------------------------------------------#

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(timestamp, " - ", ...)
}


#-----------------------------------------------------#
#-------              Import data            ---------
#-----------------------------------------------------#

log_message("Importing data...")

# Phenotype data merged with genotypes
merged_data <- readRDS(opt$data)

log_message("Defining model formula...")

# check whether covariate file provided
if (!is.null(cfile) && cfile != "" && cfile != "None") {

  covar_data <- data.table::fread(cfile)
  covariates <- colnames(covar_data[, !"IID", with = FALSE])
  covar_term <- paste(covariates, collapse = " + ")

  mformula <- paste("trait ~ geno +", covar_term) %>% as.formula()
  
  } else {
    covariates <- NULL
    mformula <- paste0("trait ~ geno") %>% as.formula()
}

log_message("Covariates: ", covariates)
log_message("Haplo.GLM model formula: ", mformula)


#-----------------------------------------------------#
#-------            Haplotype data           ---------
#-----------------------------------------------------#

log_message("Building data...")

# S1: selecting the variants
loci <- merged_data %>% dplyr::select(matches("chr[0-9]{,2}_[0-9]+_[ATCG]+_[ATCG]+"))

# list of variants
snps <- colnames(loci)

log_message(length(snps), " variants being used to build haplotypes.")

# S2: convert dosage values to major(2) or minor(1) alleles
geno_mat <- haplo.stats::geno1to2(
  round(loci, 0),
  locus.label = snps
  )

# S3: setup genotype data
geno <- haplo.stats::setupGeno(
  geno_mat,
  miss.val = c(0, NA),
  locus.label = snps
  )


#-----------------------------------------------------#
#-------       Haplotype reconstruction      ---------
#-----------------------------------------------------#

# NOTE: To have equal no. of haplotypes, NAs in phenotype must be imputed!

# Set parameters to control EM algorithm
em_ctrl <- haplo.stats::haplo.em.control(
  n.try = n_try,
  iseed = i_seed,
  insert.batch.size = n_batch,
  max.haps.limit = max_haps,
  min.posterior = min_pp  # probability of trimming off rare haplotypes
  )

#----------#
# Step 1 — Run haplotype EM once
haplo_em <- haplo.stats::haplo.em(
  geno        = geno,
  control     = em_ctrl,
  locus.label = snps
)

log_message("EM model rebuilt haplotypes using the user-specified control params.")

# Step 2 — Prepare phenotype list
phenotype_cols <- merged_data %>%
  dplyr::select(- IID, - starts_with("chr"), - all_of(covariates)) %>%
  colnames()

phenotypes <- merged_data[phenotype_cols]

log_message(length(phenotype_cols), " phenotypes being used for association test.")

#----------#
# Fitting regression model: additive haplotypes and covariate on gaussian response
# Defining Haplo.GLM model for iteration via map function
# recall parameters of haplo.GLM model from opt arguments

# Obtain haplotypes from model
grep_haplo <- function(fit) {
  summary(fit)$haplotypes %>%
    tibble::rownames_to_column(var = "Haplotype")
  }

#----------#
# Loop haplo.GLM through the phenotypes
run_haplo_glm <- function(itrait, igeno, iem, iform) {
  
  df <- data.frame(
    trait = itrait,
    haplo = igeno
  )

  # add covariates if present
  if (!is.null(covariates)) {
    df <- cbind(df, merged_data[, covariates, drop = FALSE])
  }

  # Fitting GLM model
  fit <- haplo.stats::haplo.glm(
    formula     = iform,
    family      = gaussian(),
    data        = df,
    x           = TRUE,
    haplo.em    = iem,
    locus.label = snps,
    na.action   = "na.geno.keep",
    control     = haplo.glm.control(
      haplo.freq.min = min_freq
      )
    )

    list(
      haplotype = grep_haplo(fit),
      tidy   = broom::tidy(fit) %>% suppressWarnings(),
      glance = broom::glance(fit)
  )
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

log_message("Fitting haplo.GLM...")

# Step 3 — Define a list to append model outputs
results <- vector("list", length(phenotype_cols))
names(results) <- phenotype_cols

# Step 4 — Iterating model
for (i in seq_along(phenotype_cols)) {
  
  trait_name <- phenotype_cols[i]
  trait_vec  <- merged_data[[trait_name]]
  
  results[[i]] <- run_haplo_glm(
    itrait = trait_vec,
    igeno  = geno,
    iem    = haplo_em,
    iform  = mformula
  )
  
  if (i %% 50 == 0) {
    gc()
    log_message("Processed ", i, " traits")
  }
}

log_message("Finished. Now binding list of results...")

# Step 4 — Nesting outputs in a tibble
results_tidy <- map_dfr(names(results), function(trait){
  tibble(trait_name = trait) %>%
    nest(.by = trait_name) %>% 
    dplyr::mutate(
      haplotype = list(results[[trait]]$haplotype),
      tidy      = list(results[[trait]]$tidy),
      glance    = list(results[[trait]]$glance)
    )
  }) %>%
  dplyr::select(- data)


saveRDS(results_tidy, ofile)

log_message("Stored results here: ", ofile)
