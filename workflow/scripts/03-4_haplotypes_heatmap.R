
#-----------------------------------------------------#
#-------                Logging              ---------
#-----------------------------------------------------#

# Get log path from Snakemake, fallback if missing
log_file <- tryCatch(snakemake@log[[1]], error = function(e) "logs/plot_heatmaps/default.log")

# Ensure log directory exists
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Redirect stdout & stderr to the log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # redirect stdout
sink(log_con, type = "message")  # redirect messages / stderr

# Print time & message each time
log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(timestamp, " - ", ...)
}

#-----------------------------------------------------#
#------              Load Libraries            -------
#-----------------------------------------------------#

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# Params
irds <- snakemake@input[["RDS"]]
oplt <- snakemake@output[["plt"]]
ords <- snakemake@output[["rds"]]

#-----------------------------------------------------#
#------               Import Data              -------
#-----------------------------------------------------#

log_message("Importing RDS results...")

res <- readRDS(irds)

# model coefficients
res_tidy <- res %>% 
  select(trait_name, tidy) %>%
  unnest(tidy)

#-----------------------------------------------------#
#------                 Functions              -------
#-----------------------------------------------------#

# Rename haplotypes
rename_haplo <- function(df) {
  # Keep reference haplo if identical across traits
  if (all(df$Haplotype == "Ref.")) { return("Ref.")}
  else {
    huniq <- unique(df[- c(1:2)]) #omit trait_name & Haplotype
    hash  <- apply(huniq, 1, function(x) paste(x, collapse = "_"))
    hname <- setNames(paste0("H", seq_along(hash)), hash)
    return(hname)
  }
}

#----------#
# Rename haplotype in chronological order
haplo_dict <- res %>%
  select(trait_name, haplotype) %>%
  unnest(haplotype) %>%
  select(- hap.freq) %>%
  mutate(Haplotype = str_replace_all(Haplotype, "haplo.base", "Ref.")) %>%
  filter(Haplotype != "geno.rare") %>%
  mutate(
    haplo = rename_haplo(.),
    haplo = if_else(Haplotype == "Ref.", "Ref.", haplo),
    .by = trait_name
  ) %>%
  select(trait_name, Haplotype, haplo)

#----------#
# saving tidied results for report
results_heatmap <- res_tidy %>%
  filter(str_detect(term, "geno.((\\d+)|rare)")) %>%
  # ensures to use common haplos for all traits
  inner_join(haplo_dict, join_by(trait_name, "term" == "Haplotype")) %>% 
  select(trait_name, haplo, estimate:p.value) %>%
  arrange(haplo, p.value)


# Store tidied results for HTML report
saveRDS(results_heatmap, ords)

log_message("Stored cleaned results here: ", ords)

#----------#
# convert result to wide format
results_wide <- res_tidy %>%
  dplyr::filter(str_detect(term, "geno.((\\d+)|rare)")) %>%
  # harmonize haplo names
  inner_join(haplo_dict, join_by(trait_name, "term" == "Haplotype")) %>%
  filter(haplo != "Ref.") %>%
  dplyr::select(trait_name, haplo, estimate) %>%
  pivot_wider(
    names_from = trait_name,
    values_from = estimate
    ) %>%
  # convert to matrix
  column_to_rownames("haplo") %>%
  as.matrix.data.frame()


show_dendo <- nrow(results_wide) > 1

#-----------------------------------------------------#
#------                  Heatmap               -------
#-----------------------------------------------------#

log_message("Creating heatmap...")

png(oplt, units = "in", res = 500, width = 19, height = 6)

pheatmap(
  results_wide,
  #color = hcl.colors(50, "Blue-Red 2"),
  cluster_cols = show_dendo,
  cluster_rows = show_dendo,
  clustering_method = "ward.D2",
  na_col = "white",
  border_color = NA,
  angle_col = "270"
  )

dev.off()

log_message("Stored heatmap here: ", oplt)
