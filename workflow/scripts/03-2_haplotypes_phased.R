
#-----------------------------------------------------#
# ------             Load libraries             ------
#-----------------------------------------------------#

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

#----------#
# Get arguments specified in the sbatch
option_list <- list(
  make_option("--haplo", default=NULL, help="Haplotype dosage file"),
  make_option("--pheno", default=NULL, help="Phenotype file"),
  make_option("--covariate", default=NULL, help="Covariates file (optional)"),
  make_option("--output", default=NULL, help="Output filename")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#----------#
# params
i_seed   = 777
haplo_file = opt$haplo
pheno_file = opt$pheno
covar_file = NULL #opt$covariate
ofile = opt$output


#-----------------------------------------------------#
#-------                Logging              ---------
#-----------------------------------------------------#

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(timestamp, " - ", ...)
}


# path_pipe <- "/scratch/dariush.ghasemi/projects/haplomics/"
# 
# path_haplo <- glue(path_pipe, "results/ukb_50k/phased/chr7_19k_20k_5snps.dosage.tsv")
pheno_file <- glue(path_pipe, "ukbb/pheno_50k.csv")

haplo <- read.table(haplo_file, header = TRUE, sep = "\t")
pheno <- fread(pheno_file, header = TRUE, sep = ",")

covariates <- colnames(haplo) %>% setdiff(c("sample_id", "target"))
phenotypes <- colnames(pheno[, !"IID", with = FALSE]) 

covar_term <- paste(covariates, collapse = " + ")
mformula   <- paste("trait ~", covar_term)

merged_data <- merge(pheno, haplo, by.x = "IID", by.y = "sample_id")
covar_df <- merged_data[covariates]


#----------#
# Loop haplo.GLM through the phenotypes
run_glm <- function(iform) {
  
  res <- tryCatch({
    
    # Fitting GLM model
    fit <- lm(
      iform, data = merged_data,
      #family = binomial()
    )
    
    list(
      #haplotype = grep_haplo(fit),
      tidy   = broom::tidy(fit) %>% suppressWarnings(),
      glance = broom::glance(fit)
    )
  }, warning = function(w) {
    
    message("Warning: ", conditionMessage(w))
    
    return(list(
      tidy   = NULL,
      glance = NULL,
      error  = TRUE,
      message = conditionMessage(w)
    ))
    
  }, error = function(e) {
    
    message("Error: ", conditionMessage(e))
    
    return(list(
      tidy   = NULL,
      glance = NULL,
      error  = TRUE,
      message = conditionMessage(e)
    ))
  })
  
  return(res)
}

# check whether covariate file provided
if (!is.null(covar_file) && covar_file != "" && covar_file != "None") {
  
  covar_data <- data.table::fread(covar_file)
  covariates <- colnames(covar_data[, !"IID", with = FALSE])
  covar_term <- paste(covariates, collapse = " + ")
  
  mformula <- paste("trait ~ geno +", covar_term)
  
} else {
  covariates <- NULL
  mformula <- paste0("trait ~ ", covar_term)
}


#-----------------------------------------------------#
#-------        Haplotype association        ---------
#-----------------------------------------------------#

log_message("Fitting GLM...")

# Step 3 — Define a list to append model outputs
results <- vector("list", length(phenotypes))
names(results) <- phenotypes

# Step 4 — Iterating model
for (i in seq_along(phenotypes)) {
  
  trait_name <- phenotypes[i]
  mformula <- gsub("trait", trait_name, mformula)

  results[[i]] <- run_glm(iform  = mformula)
  
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
      #haplotype = list(results[[trait]]$haplotype),
      tidy      = list(results[[trait]]$tidy),
      glance    = list(results[[trait]]$glance),
      error     = list(results[[trait]]$error),
      message   = list(results[[trait]]$message)
    )
}) %>%
  dplyr::select(- data)


saveRDS(results_tidy, ofile)

log_message("Stored results here: ", ofile)

